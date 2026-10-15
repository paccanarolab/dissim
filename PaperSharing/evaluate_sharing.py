"""Evaluates how shared publications affect the disease similarity on molecular benchmarks

Scores compared, for every pair of matrix diseases with at least one protein (UniProt mimtosp):

* original: the similarity matrix
* ablated: the matrix, with the pairs sharing publications rescored without the shared
  publications (output of ablate_shared_papers.py)
* random control: the matrix, with the same pairs rescored after removing the same number of
  random non-shared publications (mean over the repeats of ablate_shared_papers.py)
* paper-sharing baselines: shared, jaccard and adamic_adar from paper_sharing_baseline.py
  (0 for pairs sharing no publication)

Benchmarks are files of positive pairs (OMIM_ID_1 OMIM_ID_2 1); every other pair is a negative.

For each benchmark it reports ROC-AUC, average precision (AUPR) and precision at the top 100,
1,000 and 10,000 pairs, and runs a disease-level bootstrap (diseases are resampled with
replacement and each pair is weighted by the product of its two diseases' multiplicities, since
pairs sharing a disease are not independent) for:

* ablated - original and random control - original, on ROC-AUC and AUPR
* the stratified ROC-AUC test: pairs sharing publications vs pairs sharing none, on the original
  and on the ablated scores

p-values are the two-sided percentile-bootstrap ones; the smallest reportable value is 1/n_boot.
"""

import os
import sys

import numpy as np
import re
from collections import defaultdict
from sklearn.metrics import average_precision_score, roc_auc_score


def read_proteins(mimtosp):
    prots = defaultdict(set)
    with open(mimtosp, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.match(r"\s*(\d{6}):(.*)", line)
            if m:
                prots[m.group(1)].update(re.findall(r"\(([A-Z0-9]+)\)", m.group(2)))
    return prots


class SortedScores(object):
    """Scores sorted once, with ties grouped, for fast weighted ROC-AUC and AP"""

    def __init__(self, scores, iu, ju):
        self.order = np.argsort(scores, kind="stable")
        s = scores[self.order]
        self.starts = np.flatnonzero(np.r_[True, s[1:] != s[:-1]])
        self.iu, self.ju = iu[self.order], ju[self.order]

    def sorted(self, values):
        return values[self.order]

    def weights(self, counts):
        return counts[self.iu] * counts[self.ju]

    def auc(self, y_o, w_o):
        wp = np.add.reduceat(w_o * y_o, self.starts)
        wn = np.add.reduceat(w_o * ~y_o, self.starts)
        below = np.cumsum(wn) - wn
        return float((wp * (below + 0.5 * wn)).sum() / (wp.sum() * wn.sum()))

    def ap(self, y_o, w_o):
        tp = np.add.reduceat(w_o * y_o, self.starts)[::-1]
        fp = np.add.reduceat(w_o * ~y_o, self.starts)[::-1]
        ctp, cfp = np.cumsum(tp), np.cumsum(fp)
        with np.errstate(invalid="ignore", divide="ignore"):
            precision = np.where(ctp + cfp > 0, ctp / (ctp + cfp), 0.0)
        return float((tp * precision).sum() / tp.sum())


def precision_at(scores, y, k):
    top = np.argpartition(-scores, k)[:k]
    return float(y[top].mean())


def bootstrap_summary(values, point):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    p = min(1.0, 2 * min((values <= 0).mean(), (values >= 0).mean()))
    return point, np.percentile(values, 2.5), np.percentile(values, 97.5), max(p, 1 / len(values)), p == 0


def run(args):
    ids = [line.strip() for line in open(args.matrix_ids) if line.strip()]
    prots = read_proteins(args.mimtosp)
    keep = [i for i, omim in enumerate(ids) if prots.get(omim)]
    universe = [ids[i] for i in keep]
    n = len(universe)
    position = {omim: k for k, omim in enumerate(universe)}
    print(f"{len(ids)} matrix diseases, {n} with at least one protein -> {n * (n - 1) // 2:,} pairs")

    matrix = np.load(args.matrix, mmap_mode="r")
    iu, ju = np.triu_indices(n, 1)
    iu, ju = iu.astype(np.int32), ju.astype(np.int32)
    original = np.asarray(matrix[np.ix_(keep, keep)])[iu, ju].astype(np.float64)
    num_pairs = len(original)

    def pair_index(a, b):
        i, j = sorted((position[a], position[b]))
        return i * n - i * (i + 1) // 2 + (j - i - 1)

    def read_pairs(filename, columns):
        rows = []
        with open(filename) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if fields[0] in position and fields[1] in position and fields[0] != fields[1]:
                    rows.append([pair_index(fields[0], fields[1])] + [float(fields[c]) for c in columns])
        return np.array(rows)

    # ablation: columns original(5) ablated(6) ablated_empty(7) random_mean(8)
    ablation = read_pairs(args.ablation, [5, 6, 7, 8])
    idx = ablation[:, 0].astype(np.int64)
    largest = np.abs(original[idx] - ablation[:, 1]).max()
    if largest > 1e-6:
        sys.exit(f"The ablation file was computed from a different matrix (largest difference {largest:.3g})")
    shares = np.zeros(num_pairs, bool)
    shares[idx] = True
    ablated = original.copy()
    ablated[idx] = ablation[:, 2]
    random_control = original.copy()
    random_control[idx] = ablation[:, 4]
    print(f"Pairs sharing publications: {len(idx):,} ({int(ablation[:, 3].sum())} with no publication left "
          f"after removing the shared ones)")

    baseline = read_pairs(args.baseline, [2, 3, 4])
    baselines = {}
    for k, name in enumerate(["shared", "jaccard", "adamic_adar"]):
        values = np.zeros(num_pairs)
        values[baseline[:, 0].astype(np.int64)] = baseline[:, k + 1]
        baselines[f"baseline: {name}"] = values

    labels = {}
    for benchmark in args.benchmarks:
        y = np.zeros(num_pairs, bool)
        with open(benchmark) as f:
            for line in f:
                fields = line.split()
                if len(fields) >= 2 and fields[0] in position and fields[1] in position and fields[0] != fields[1]:
                    y[pair_index(fields[0], fields[1])] = True
        labels[os.path.basename(benchmark)] = y

    scores = {"original": original, "ablated": ablated, "random control": random_control, **baselines}

    # point estimates
    rows = []
    for bench, y in labels.items():
        for name, s in scores.items():
            rows.append((bench, name, "all pairs", roc_auc_score(y, s), average_precision_score(y, s), y.mean(),
                         int(y.sum()), *(precision_at(s, y, k) for k in (100, 1000, 10000))))
        for name in ("original", "ablated"):
            for stratum, mask in (("shares publications", shares), ("shares none", ~shares)):
                s = scores[name][mask]
                rows.append((bench, name, stratum, roc_auc_score(y[mask], s), average_precision_score(y[mask], s),
                             y[mask].mean(), int(y[mask].sum()), np.nan, np.nan, np.nan))

    # disease-level bootstrap
    print(f"Bootstrap ({args.n_boot} resamples)..")
    sorted_scores = {name: SortedScores(scores[name], iu, ju) for name in ("original", "ablated", "random control")}
    sorted_labels = {(name, bench): ss.sorted(y) for name, ss in sorted_scores.items() for bench, y in labels.items()}
    sorted_shares = {name: ss.sorted(shares) for name, ss in sorted_scores.items()}
    rng = np.random.default_rng(args.seed)
    boot = defaultdict(list)
    for r in range(args.n_boot):
        counts = np.bincount(rng.integers(0, n, n), minlength=n).astype(np.float64)
        for name, ss in sorted_scores.items():
            w = ss.weights(counts)
            for bench in labels:
                y_o = sorted_labels[(name, bench)]
                boot[(bench, name, "auc")].append(ss.auc(y_o, w))
                boot[(bench, name, "ap")].append(ss.ap(y_o, w))
                if name in ("original", "ablated"):
                    with np.errstate(invalid="ignore", divide="ignore"):
                        boot[(bench, name, "auc share")].append(ss.auc(y_o, w * sorted_shares[name]))
                        boot[(bench, name, "auc no share")].append(ss.auc(y_o, w * ~sorted_shares[name]))
        if (r + 1) % 20 == 0:
            print(f"  {r + 1}/{args.n_boot}", flush=True)

    point = {(r[0], r[1], r[2]): r for r in rows}
    tests = []
    for bench in labels:
        for name in ("ablated", "random control"):
            for metric, col in (("ROC-AUC", 3), ("AUPR", 4)):
                key = "auc" if metric == "ROC-AUC" else "ap"
                diff = np.array(boot[(bench, name, key)]) - np.array(boot[(bench, "original", key)])
                tests.append((bench, f"{name} - original", metric,
                              *bootstrap_summary(diff, point[(bench, name, "all pairs")][col]
                                                 - point[(bench, "original", "all pairs")][col])))
        for name in ("original", "ablated"):
            diff = np.array(boot[(bench, name, "auc share")]) - np.array(boot[(bench, name, "auc no share")])
            tests.append((bench, f"{name}: shares publications - shares none", "ROC-AUC",
                          *bootstrap_summary(diff, point[(bench, name, "shares publications")][3]
                                             - point[(bench, name, "shares none")][3])))

    with open(args.outfile, "w") as o:
        o.write("#benchmark\tscore\tpairs\troc_auc\taupr\tbase_rate\tpositives\tp@100\tp@1000\tp@10000\n")
        for r in rows:
            o.write("\t".join(str(v) for v in r) + "\n")
        o.write("#benchmark\ttest\tmetric\tdifference\tci_low\tci_high\tp_value\tp_below_resolution\n")
        for t in tests:
            o.write("\t".join(str(v) for v in t) + "\n")

    for bench in labels:
        print(f"\n===== {bench}: {int(labels[bench].sum()):,} positives, base rate {100 * labels[bench].mean():.3f}%")
        print(f"  {'score':22s} {'pairs':20s} {'ROC-AUC':>8s} {'AUPR':>7s} {'base%':>7s} {'#pos':>7s} "
              f"{'P@100':>6s} {'P@1k':>6s} {'P@10k':>6s}")
        for r in rows:
            if r[0] == bench:
                pk = "" if np.isnan(r[7]) else f" {r[7]:6.2f} {r[8]:6.3f} {r[9]:6.3f}"
                print(f"  {r[1]:22s} {r[2]:20s} {r[3]:8.4f} {r[4]:7.4f} {100 * r[5]:7.3f} {r[6]:7,}{pk}")
        print(f"  bootstrap tests ({args.n_boot} disease resamples):")
        for t in tests:
            if t[0] == bench:
                print(f"    {t[1]:45s} {t[2]:8s} {t[3]:+.4f}  95% CI {t[4]:+.4f} .. {t[5]:+.4f}  "
                      f"p {'<' if t[7] else '='} {t[6]:.3f}")


if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    req = aparser.add_argument_group("required arguments")
    req.add_argument("--matrix", help="dense similarity matrix (.npy)", required=True)
    req.add_argument("--matrix-ids", help="OMIM ids of the matrix rows, one per line", required=True)
    req.add_argument("--ablation", help="output of ablate_shared_papers.py", required=True)
    req.add_argument("--baseline", help="output of paper_sharing_baseline.py", required=True)
    req.add_argument("--mimtosp", help="UniProt mimtosp.txt (OMIM to Swiss-Prot)", required=True)
    req.add_argument("--benchmarks", nargs="+", help="benchmark files of positive pairs", required=True)
    req.add_argument("--outfile", "--o", help="path to the output results file", required=True)
    aparser.add_argument("--n-boot", type=int, default=200, help="bootstrap resamples")
    aparser.add_argument("--seed", type=int, default=0)
    run(aparser.parse_args())
