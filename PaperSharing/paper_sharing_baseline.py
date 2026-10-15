"""Cheap paper-sharing baseline between OMIM entries

Scores every pair of OMIM entries by the publications they both reference, using only the
omim2pubmed mapping (no MeSH, no ontology). Pairs that share no publication score 0 and are
not written, so the output is sparse: any pair missing from it has every score equal to 0.

Scores written for each pair sharing at least one publication:

* shared: number of shared publications
* jaccard: shared / size of the union of both reference lists
* adamic_adar: sum over shared publications p of 1 / log(1 + deg(p)), where deg(p) is the number
  of OMIM entries referencing p. A publication referenced by few entries counts more, the same
  idea as the information content used by the MeSH measures.

Output format (tab separated, with a header line starting with '#'):

    OMIM_ID_1  OMIM_ID_2  shared  jaccard  adamic_adar
"""

import itertools
import math
from collections import defaultdict

import rich.progress


def read_omim2pubmed(filename):
    papers = {}
    with open(filename, "r") as f:
        for line in f:
            fields = line.split()
            if len(fields) > 1:
                papers[fields[0]] = set(fields[1:])
    return papers


def compute(papers):
    cited_by = defaultdict(list)
    for omim, pubmeds in papers.items():
        for pubmed in pubmeds:
            cited_by[pubmed].append(omim)
    shared = defaultdict(int)
    adamic_adar = defaultdict(float)
    for pubmed, omims in rich.progress.track(cited_by.items(), description="Finding shared publications..."):
        if len(omims) < 2:
            continue
        weight = 1.0 / math.log(1 + len(omims))
        for pair in itertools.combinations(sorted(omims), 2):
            shared[pair] += 1
            adamic_adar[pair] += weight
    for (a, b), n in shared.items():
        yield a, b, n, n / len(papers[a] | papers[b]), adamic_adar[(a, b)]


def run(omim2pubmed, outfile):
    papers = read_omim2pubmed(omim2pubmed)
    with open(outfile, "w") as o:
        o.write("#OMIM_ID_1\tOMIM_ID_2\tshared\tjaccard\tadamic_adar\n")
        for a, b, n, jaccard, adamic_adar in compute(papers):
            o.write(f"{a}\t{b}\t{n}\t{jaccard:.6g}\t{adamic_adar:.6g}\n")


if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    req = aparser.add_argument_group("required arguments")
    req.add_argument("--omim-to-pubmed", "--in", help="OMIM2PubMed mapping file", required=True)
    req.add_argument("--outfile", "--o", help="path to the output file", required=True)
    args = aparser.parse_args()
    run(args.omim_to_pubmed, args.outfile)
