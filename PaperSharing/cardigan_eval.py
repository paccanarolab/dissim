"""Compares disease similarities through Cardigan disease-gene prediction

For every disease with known genes (OMIM morbidmap) we predict its genes the way Cardigan does for
a molecularly uncharacterised disease: no seed genes are given, the only input is the disease, so
the initial labels come entirely from the genes of similar diseases (weighted by Cardigan's
sigmoid of the similarity) and are then diffused on the PPI. The known genes of the disease are
the targets. This is the setting where the result depends entirely on the similarity.

Similarities compared (all restricted to the same set of diseases, so they draw from the same
pool of similar diseases):

* caniza_2017: the similarity shipped with Cardigan
* caniza_new: the new similarity matrix
* caniza_new_ablated: caniza_new, with the pairs sharing publications rescored without the shared
  publications (output of ablate_shared_papers.py)
* paper_sharing: the adamic_adar score of paper_sharing_baseline.py, rank-mapped onto the values
  of caniza_new (its K non-zero pairs get the K largest caniza_new values, in the same order), so
  that Cardigan's sigmoid, tuned on the Caniza scale, treats it comparably

Everything else (PPI loading, morbidmap, sigmoid, diffusion) is Cardigan's code, unmodified.

Per disease and similarity, from the diffusion scores of all PPI genes:

* auroc: ROC-AUC of the targets against every other gene (ties averaged)
* recall@k: fraction of targets ranked in the top k genes
* mean_rank: mean rank of the targets (1 is best)

A disease with no similar disease (nothing to diffuse) is scored as a random ranking (auroc 0.5,
recall 0, mean_rank (n + 1) / 2) and flagged.

Cardigan's ConsistencyMethod.diffuse is replaced by the same linear solve without its
`if self.S == []` guard, which recent SciPy versions reject for sparse matrices.
"""

import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.stats import rankdata, wilcoxon

RECALL_AT = (10, 50, 100, 200)
_state = {}


class MatrixSimilarity(object):
    """Drop-in replacement for Cardigan's DISSIMUtilities.DiseaseSimilarity backed by a .npy matrix"""

    def __init__(self, matrix, ids, diseases, overrides=None):
        self.matrix = np.load(matrix, mmap_mode="r")
        self.index = {m: i for i, m in enumerate(ids)}
        self.diseases = set(diseases)
        self.overrides = overrides or {}

    def getDiseases(self):
        return self.diseases

    def similarity(self, mimA, mimB):
        pair = (min(mimA, mimB), max(mimA, mimB))
        if pair in self.overrides:
            return self.overrides[pair]
        return float(self.matrix[self.index[mimA], self.index[mimB]])


class SparseSimilarity(object):
    """Similarity given only for some pairs, 0 elsewhere"""

    def __init__(self, pairs, diseases):
        self.pairs = pairs
        self.diseases = set(diseases)

    def getDiseases(self):
        return self.diseases

    def similarity(self, mimA, mimB):
        return self.pairs.get((min(mimA, mimB), max(mimA, mimB)), 0)


def read_ids(filename):
    return [line.strip() for line in open(filename) if line.strip()]


def read_pair_values(filename, column):
    values = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            values[(min(fields[0], fields[1]), max(fields[0], fields[1]))] = float(fields[column])
    return values


def rank_map_baseline(args, universe):
    """adamic_adar values -> the largest caniza_new values, preserving the order (ties share a value)"""
    baseline = {p: v for p, v in read_pair_values(args.baseline, 4).items()
                if p[0] in universe and p[1] in universe and p[0] != p[1]}
    ids = read_ids(args.new_ids)
    keep = sorted(i for i, m in enumerate(ids) if m in universe)
    matrix = np.load(args.new_matrix, mmap_mode="r")
    sub = np.asarray(matrix[np.ix_(keep, keep)])
    upper = sub[np.triu_indices(len(keep), 1)]
    k = len(baseline)
    target = np.sort(np.partition(upper, len(upper) - k)[len(upper) - k:])
    pairs = list(baseline)
    ranks = rankdata([baseline[p] for p in pairs], method="average") - 1
    mapped = np.interp(ranks, np.arange(k), target)
    return dict(zip(pairs, mapped.tolist()))


def init_worker(args, universe, baseline_mapped):
    sys.path.insert(0, os.path.join(args.cardigan, "lib"))
    import OMIMUtilities
    import PPIUtilities
    import diffusionUtilities
    import scipy.sparse.linalg as la

    def consistency_diffuse(self, initialGuess=[]):
        # ConsistencyMethod.diffuse without its `if self.S == []` guard, which recent SciPy rejects
        # for sparse matrices; the linear system solved is the same
        if len(initialGuess) > 0:
            self.Y = initialGuess
        return la.cg(self.IalphaS, self.Y)[0]

    diffusionUtilities.ConsistencyMethod.diffuse = consistency_diffuse

    ppi = PPIUtilities.PPI(args.ppi, columns="hippie", ppiName="HIPPIE")
    mim2gene = OMIMUtilities.Mim2Gene(args.mim2gene)
    morbidmap = OMIMUtilities.MorbidMap(args.morbidmap, separator="\t", mim2gene=mim2gene, verbose=False)

    ablated = {(min(a, b), max(a, b)): v for (a, b), v in read_pair_values(args.ablation, 6).items()}
    similarities = {
        "caniza_2017": MatrixSimilarity(args.old_matrix, read_ids(args.old_ids), universe),
        "caniza_new": MatrixSimilarity(args.new_matrix, read_ids(args.new_ids), universe),
        "caniza_new_ablated": MatrixSimilarity(args.new_matrix, read_ids(args.new_ids), universe, ablated),
        "paper_sharing": SparseSimilarity(baseline_mapped, universe),
    }
    # same parameters Cardigan's Config.setPPI('hippie') sets
    diffusion = diffusionUtilities.Diffusion(ppi, morbidmap, similarities["caniza_new"], verbose=False)
    diffusion.setSigmoidScale(3.0, 4.0, 0.6)
    diffusion.setDiffusionParameters("consistency", alpha=0.75)
    _state.update(diffusion=diffusion, similarities=similarities, ppi=ppi, morbidmap=morbidmap,
                  gene_index={g: i for i, g in enumerate(ppi.vertices)})


def evaluate_disease(mim):
    diffusion, ppi = _state["diffusion"], _state["ppi"]
    targets = sorted(set(_state["morbidmap"].getGenes(mim)) & set(ppi.vertices))
    target_idx = np.array([_state["gene_index"][g] for g in targets])
    out = []
    for name, similarity in _state["similarities"].items():
        diffusion.ds = similarity
        diffusion.setDiffusionMims([mim])
        diffusion.setDiffusionGenes([])
        diffusion.setDiffusionTargets(targets)
        result = diffusion.diffuse()
        n_similar = sum(1 for w in (diffusion.lastGeneWeights or {}).values() if w > 0)
        if not result:
            # no input: score as a random ranking (expected AUROC 0.5, expected rank (n + 1) / 2)
            out.append((mim, name, len(targets), 0.5, *[0.0] * len(RECALL_AT), (len(ppi.vertices) + 1) / 2, 0, 1))
            continue
        scores = np.asarray(diffusion.lastDiffusionVector).ravel()
        # rank 1 = highest score, ties get the average rank
        ranks = rankdata(-scores, method="average")
        t_ranks = ranks[target_idx]
        n, n_t = len(scores), len(target_idx)
        # AUC = P(target scored above a non-target), from the ranks among all genes
        auc = 1 - (t_ranks.sum() - n_t * (n_t + 1) / 2) / (n_t * (n - n_t))
        out.append((mim, name, n_t, auc, *[(t_ranks <= k).mean() for k in RECALL_AT], t_ranks.mean(),
                    n_similar, 0))
    return out


def run(args):
    t0 = time.time()
    sys.path.insert(0, os.path.join(args.cardigan, "lib"))
    import OMIMUtilities
    import PPIUtilities

    # diseases available to every similarity
    universe = set(read_ids(args.old_ids)) & set(read_ids(args.new_ids))
    ppi = PPIUtilities.PPI(args.ppi, columns="hippie", ppiName="HIPPIE")
    mim2gene = OMIMUtilities.Mim2Gene(args.mim2gene)
    morbidmap = OMIMUtilities.MorbidMap(args.morbidmap, separator="\t", mim2gene=mim2gene, verbose=False)
    vertices = set(ppi.vertices)
    diseases = sorted(m for m in set(morbidmap.getDiseases()) & universe if set(morbidmap.getGenes(m)) & vertices)
    if args.limit:
        diseases = diseases[:args.limit]
    print(f"PPI: {len(vertices)} genes | diseases in both similarities: {len(universe)} | "
          f"evaluated (in morbidmap, with genes in the PPI): {len(diseases)}", flush=True)

    baseline_mapped = rank_map_baseline(args, universe)
    print(f"paper_sharing: {len(baseline_mapped)} non-zero pairs, mapped onto caniza_new values "
          f"{min(baseline_mapped.values()):.3f}..{max(baseline_mapped.values()):.3f}", flush=True)

    rows = []
    with Pool(args.workers, initializer=init_worker, initargs=(args, universe, baseline_mapped)) as pool:
        for k, result in enumerate(pool.imap_unordered(evaluate_disease, diseases, chunksize=4)):
            rows.extend(result)
            if (k + 1) % 200 == 0:
                print(f"  {k + 1}/{len(diseases)} diseases, {time.time() - t0:.0f}s", flush=True)

    header = ["mim", "similarity", "targets", "auroc", *[f"recall@{k}" for k in RECALL_AT], "mean_rank",
              "similar_genes", "nothing_to_diffuse"]
    with open(args.outfile, "w") as o:
        o.write("\t".join(header) + "\n")
        for r in sorted(rows):
            o.write("\t".join(str(v) for v in r) + "\n")

    names = ["caniza_2017", "caniza_new", "caniza_new_ablated", "paper_sharing"]
    table = {name: {r[0]: r for r in rows if r[1] == name} for name in names}
    print(f"\n{len(diseases)} diseases, {time.time() - t0:.0f}s")
    print(f"{'similarity':20s} {'AUROC':>7s} " + " ".join(f"{'R@' + str(k):>6s}" for k in RECALL_AT)
          + f" {'mean rank':>9s} {'no input':>8s}")
    for name in names:
        r = np.array([table[name][m][3:3 + 1 + len(RECALL_AT)] for m in diseases], dtype=float)
        mean_rank = np.nanmean([table[name][m][3 + 1 + len(RECALL_AT)] for m in diseases])
        empty = sum(table[name][m][-1] for m in diseases)
        print(f"{name:20s} {r[:, 0].mean():7.4f} " + " ".join(f"{v:6.3f}" for v in r[:, 1:].mean(0))
              + f" {mean_rank:9.1f} {empty:8d}")
    print("\npaired Wilcoxon signed-rank tests on per-disease AUROC:")
    for a, b in [("caniza_new", "caniza_2017"), ("caniza_new_ablated", "caniza_new"),
                 ("paper_sharing", "caniza_new"), ("paper_sharing", "caniza_2017")]:
        x = np.array([table[a][m][3] for m in diseases]); y = np.array([table[b][m][3] for m in diseases])
        d = x - y
        p = wilcoxon(x, y).pvalue if np.any(d != 0) else 1.0
        print(f"  {a} - {b}: mean {d.mean():+.4f}, better/worse/equal {int((d > 0).sum())}/{int((d < 0).sum())}"
              f"/{int((d == 0).sum())}, p = {p:.2g}")


if __name__ == "__main__":
    import argparse

    DISSIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    DESKTOP = os.path.join(DISSIM, "..")
    aparser = argparse.ArgumentParser()
    aparser.add_argument("--cardigan", default=os.path.join(DESKTOP, "Cardigan"))
    aparser.add_argument("--ppi", default=os.path.join(DESKTOP, "RareDisWP1", "BuildBenchmark", "files", "hippie_current.txt"))
    aparser.add_argument("--morbidmap", default=os.path.join(DESKTOP, "Cardigan", "data", "mappings", "omim", "morbidmap.txt"))
    aparser.add_argument("--mim2gene", default=os.path.join(DESKTOP, "Cardigan", "data", "mappings", "omim", "mim2gene.txt"))
    aparser.add_argument("--old-matrix", default=os.path.join(DISSIM, "data", "cardigan", "caniza_2017.npy"))
    aparser.add_argument("--old-ids", default=os.path.join(DISSIM, "data", "cardigan", "caniza_2017_ids.txt"))
    aparser.add_argument("--new-matrix", default=os.path.join(DISSIM, "data", "similarity_run", "resnik_two.npy"))
    aparser.add_argument("--new-ids", default=os.path.join(DISSIM, "data", "similarity_run", "resnik_two_ids.txt"))
    aparser.add_argument("--ablation", default=os.path.join(DISSIM, "data", "ablation_resnik_two.tsv"))
    aparser.add_argument("--baseline", default=os.path.join(DISSIM, "data", "paper_sharing.tsv"))
    aparser.add_argument("--outfile", "--o", default=os.path.join(DISSIM, "data", "cardigan", "cardigan_eval.tsv"))
    aparser.add_argument("--workers", type=int, default=6)
    aparser.add_argument("--limit", type=int, default=0, help="only evaluate the first N diseases (testing)")
    run(aparser.parse_args())
