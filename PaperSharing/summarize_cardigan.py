"""Summarises cardigan_eval.py results into small CSV files (summary, paired tests, per disease, run)

A disease "shares papers" when at least one other disease of the evaluated similarity universe
references a publication it references too: only those diseases can change when the shared
publications are removed (caniza_new_ablated) or when only paper sharing is used.
"""

import csv
import os

import numpy as np
from scipy.stats import wilcoxon

NAMES = ["caniza_2017", "caniza_new", "caniza_new_ablated", "paper_sharing"]
LABELS = {"caniza_2017": "Caniza 2017", "caniza_new": "Caniza new",
          "caniza_new_ablated": "Caniza new, shared papers removed", "paper_sharing": "Paper sharing"}
METRICS = ["auroc", "recall@10", "recall@50", "recall@100", "recall@200", "mean_rank"]
COMPARISONS = [("caniza_new", "caniza_2017"), ("caniza_new_ablated", "caniza_new"),
               ("paper_sharing", "caniza_new"), ("paper_sharing", "caniza_2017")]


def run(args):
    rows = list(csv.DictReader(open(args.results), delimiter="\t"))
    data = {}
    for r in rows:
        data.setdefault(r["mim"], {})[r["similarity"]] = r
    mims = sorted(m for m in data if all(n in data[m] for n in NAMES))

    universe = {l.strip() for l in open(args.old_ids)} & {l.strip() for l in open(args.new_ids)}
    sharing = set()
    with open(args.baseline) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a, b = line.split()[:2]
            if a in universe and b in universe:
                sharing.update((a, b))
    shares = np.array([m in sharing for m in mims])

    def values(name, metric):
        return np.array([float(data[m][name][metric]) for m in mims])

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "summary.csv"), "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["similarity", "label", "diseases", "auroc_mean", "auroc_median", "recall_10", "recall_50",
                    "recall_100", "recall_200", "mean_rank", "no_input"])
        for n in NAMES:
            a = values(n, "auroc")
            w.writerow([n, LABELS[n], len(mims), round(a.mean(), 5), round(float(np.median(a)), 5),
                        *[round(values(n, f"recall@{k}").mean(), 5) for k in (10, 50, 100, 200)],
                        round(values(n, "mean_rank").mean(), 1),
                        int(values(n, "nothing_to_diffuse").sum())])

    rng = np.random.default_rng(0)
    with open(os.path.join(args.outdir, "tests.csv"), "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["comparison", "a", "b", "diseases_subset", "diseases", "mean_diff", "ci_low", "ci_high",
                    "better", "worse", "equal", "p_value"])
        for a, b in COMPARISONS:
            for subset, mask in (("all", np.ones(len(mims), bool)), ("shares_papers", shares)):
                d = (values(a, "auroc") - values(b, "auroc"))[mask]
                boots = [rng.choice(d, len(d)).mean() for _ in range(2000)]
                p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
                w.writerow([f"{LABELS[a]} vs {LABELS[b]}", a, b, subset, len(d), f"{d.mean():.5f}",
                            f"{np.percentile(boots, 2.5):.5f}", f"{np.percentile(boots, 97.5):.5f}",
                            int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum()), f"{p:.3g}"])

    with open(os.path.join(args.outdir, "per_disease.csv"), "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["mim", "targets", "shares_papers"] + [f"auroc_{n}" for n in NAMES]
                   + [f"recall100_{n}" for n in NAMES])
        for m, s in zip(mims, shares):
            w.writerow([m, data[m][NAMES[0]]["targets"], "yes" if s else "no"]
                       + [f"{float(data[m][n]['auroc']):.5f}" for n in NAMES]
                       + [f"{float(data[m][n]['recall@100']):.3f}" for n in NAMES])

    with open(os.path.join(args.outdir, "run.csv"), "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["key", "setting", "value"])
        for k, label, v in args.run_info:
            w.writerow([k, label, v])
        w.writerow(["diseases", "Diseases evaluated", len(mims)])
        w.writerow(["diseases_sharing", "Diseases that share papers", int(shares.sum())])
    print(f"{len(mims)} diseases, {int(shares.sum())} share papers -> {args.outdir}")


if __name__ == "__main__":
    import argparse

    DISSIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    aparser = argparse.ArgumentParser()
    aparser.add_argument("--results", default=os.path.join(DISSIM, "data", "cardigan", "cardigan_eval.tsv"))
    aparser.add_argument("--old-ids", default=os.path.join(DISSIM, "data", "cardigan", "caniza_2017_ids.txt"))
    aparser.add_argument("--new-ids", default=os.path.join(DISSIM, "data", "similarity_run", "resnik_two_ids.txt"))
    aparser.add_argument("--baseline", default=os.path.join(DISSIM, "data", "paper_sharing.tsv"))
    aparser.add_argument("--log", default=os.path.join(DISSIM, "data", "cardigan", "cardigan_eval.log"))
    aparser.add_argument("--outdir", default=os.path.join(DISSIM, "data", "cardigan", "dashboard"))
    args = aparser.parse_args()
    # read from the cardigan_eval.py log
    log = open(args.log).read()
    import re
    genes, universe_size, _ = re.search(r"PPI: (\d+) genes \| diseases in both similarities: (\d+) \| evaluated.*?: (\d+)", log).groups()
    pairs = re.search(r"paper_sharing: (\d+) non-zero pairs", log).group(1)
    args.run_info = [("ppi", "PPI network", "HIPPIE (current)"), ("ppi_genes", "PPI genes", genes),
                     ("universe", "Diseases in both Caniza matrices", universe_size),
                     ("sharing_pairs", "Paper-sharing pairs used", pairs),
                     ("alpha", "Diffusion alpha", 0.75), ("sig10", "Sigmoid 10% point", 3.0),
                     ("sig90", "Sigmoid 90% point", 4.0), ("sigh", "Sigmoid height", 0.6)]
    run(args)
