"""Recomputes the Resnik (MAX) disease similarity of every pair of OMIM entries that share
publications, after removing the shared publications from both entries.

With the term-wise MAX strategy used by compute_combined_similarity.py, the similarity of two
diseases is the largest Resnik similarity between any of their MeSH descriptors, and the Resnik
similarity of two descriptors is the largest information content (IC) among their common
ancestors (each descriptor is its own ancestor). Both maxima combine into one:

    sim(A, B) = max IC(c) over the descriptors c annotating both A and B (after up-propagation)

so a single pair can be recomputed from its two annotation sets, without the descriptor matrix.

The IC values are always those of the full corpus, so the only thing that changes for a pair is
which publications its two diseases contribute. A pair where one disease has no publication left
gets a similarity of 0.

As a control, the same number of publications is removed at random from the publications each
disease does NOT share, keeping the shared ones (--random-repeats times). If a disease has fewer
non-shared publications than needed, all of them are removed and the pair is flagged.

Before writing anything, the full (non-ablated) scores are recomputed and compared with the
similarity matrix, to check that the result reproduces compute_combined_similarity.py.
"""

import itertools
import math
import os
import random
import sys
from collections import defaultdict

import numpy as np
import rich.progress

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ComputeSimilarities"))
from mesh_parser import MeSHParser  # noqa: E402
from annotation import AnnotationParser  # noqa: E402

CATEGORIES = {
    "TWO": ["A", "C"],
    "FIVE": ["A", "C", "D", "E", "G"],
    "ALL": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "V", "Z"],
}


def read_mapping(filename):
    # same semantics as MIM2MESH.readMappingFile: one line per key, the last one wins
    values = {}
    with open(filename, "r") as f:
        for line in f:
            fields = line.split()
            if fields:
                values[fields[0]] = fields[1:]
    return values


class Scorer(object):

    def __init__(self, mesh_file, categories, omim2mesh, pubmed2mesh):
        print("Loading MeSH and the annotation (same code as compute_combined_similarity.py)..")
        thesaurus = MeSHParser(mesh_file, CATEGORIES[categories]).get_thesaurus(True)
        annotation = AnnotationParser(thesaurus, omim2mesh).get_annotations()
        self.annotation = annotation
        num_objects = len(annotation.get_objects())
        self.ic = {d: -math.log10(annotation.num_annot_per_descriptor(d) / num_objects)
                   for d in annotation.get_descriptors()}
        self.in_thesaurus = set(thesaurus.get_node_ids())
        self.ancestors = {}
        self.thesaurus = thesaurus
        self.pubmed2mesh = pubmed2mesh
        self.paper_cache = {}
        print(f"{num_objects} annotated diseases, {len(self.ic)} annotated descriptors")

    def __ancestors(self, descriptor):
        if descriptor not in self.ancestors:
            self.ancestors[descriptor] = frozenset(n.get_identifier()
                                                   for n in self.thesaurus.get_node(descriptor).get_ancestors())
        return self.ancestors[descriptor]

    def paper_terms(self, pubmed):
        """Up-propagated descriptors contributed by one publication"""
        if pubmed not in self.paper_cache:
            terms = set()
            for descriptor in set(self.pubmed2mesh.get(pubmed, [])) & self.in_thesaurus:
                terms |= self.__ancestors(descriptor)
            self.paper_cache[pubmed] = frozenset(terms)
        return self.paper_cache[pubmed]

    def disease_terms(self, pubmeds):
        terms = set()
        for pubmed in pubmeds:
            terms |= self.paper_terms(pubmed)
        return terms

    def similarity(self, terms_a, terms_b):
        common = terms_a & terms_b
        return max((self.ic[c] for c in common), default=0.0)


def run(args):
    omim2pubmed = {k: set(v) for k, v in read_mapping(args.omim_to_pubmed).items()}
    pubmed2mesh = read_mapping(args.pubmed_to_mesh)
    scorer = Scorer(args.mesh, args.categories, args.omim_to_mesh, pubmed2mesh)

    ids = [line.strip() for line in open(args.matrix_ids) if line.strip()]
    index = {omim: i for i, omim in enumerate(ids)}
    matrix = np.load(args.matrix, mmap_mode="r")

    # 1. the annotation rebuilt from the publications must be the one used by the pipeline
    mismatched = [omim for omim in ids
                  if scorer.disease_terms(omim2pubmed.get(omim, ())) != scorer.annotation.get_descriptors_per_object(omim)]
    print(f"Diseases whose annotation differs when rebuilt from their publications: {len(mismatched)}")

    # 2. pairs of matrix diseases sharing publications
    cited_by = defaultdict(list)
    for omim in ids:
        for pubmed in omim2pubmed.get(omim, ()):
            cited_by[pubmed].append(omim)
    pairs = set()
    for omims in cited_by.values():
        pairs.update(itertools.combinations(sorted(omims), 2))
    pairs = sorted(pairs)
    print(f"Pairs of matrix diseases sharing at least one publication: {len(pairs)}")

    # 3. check that the full scores reproduce the matrix (sharing pairs + a random sample of other pairs)
    rng = random.Random(args.seed)
    sample = pairs + [tuple(rng.sample(ids, 2)) for _ in range(args.validation_sample)]
    full_terms = {omim: scorer.annotation.get_descriptors_per_object(omim) for omim in ids}
    diffs = np.array([abs(scorer.similarity(full_terms[a], full_terms[b]) - matrix[index[a], index[b]])
                      for a, b in rich.progress.track(sample, description="Validating against the matrix...")])
    print(f"Validation on {len(sample)} pairs: max |recomputed - matrix| = {diffs.max():.3g}, "
          f"pairs differing by more than 1e-9: {(diffs > 1e-9).sum()}")
    if diffs.max() > args.tolerance:
        sys.exit("The recomputed scores do not reproduce the matrix: check the MeSH file, the categories "
                 "and that omim2mesh/omim2pubmed/pubmed2mesh are the ones the matrix was computed from.")

    # 4. ablation and random control
    with open(args.outfile, "w") as o:
        o.write("#OMIM_ID_1\tOMIM_ID_2\tshared\tpapers_1\tpapers_2\toriginal\tablated\tablated_empty"
                f"\trandom_mean\trandom_sd\trandom_capped\n")
        for a, b in rich.progress.track(pairs, description="Removing shared publications..."):
            papers_a, papers_b = omim2pubmed[a], omim2pubmed[b]
            shared = papers_a & papers_b
            only_a, only_b = papers_a - shared, papers_b - shared
            ablated = scorer.similarity(scorer.disease_terms(only_a), scorer.disease_terms(only_b))
            empty = int(not only_a or not only_b)

            k = len(shared)
            capped = int(len(only_a) < k or len(only_b) < k)
            pair_rng = random.Random(f"{args.seed}-{a}-{b}")
            list_a, list_b = sorted(only_a), sorted(only_b)
            random_scores = []
            for _ in range(args.random_repeats):
                drop_a = set(pair_rng.sample(list_a, min(k, len(list_a))))
                drop_b = set(pair_rng.sample(list_b, min(k, len(list_b))))
                random_scores.append(scorer.similarity(scorer.disease_terms(papers_a - drop_a),
                                                       scorer.disease_terms(papers_b - drop_b)))
            o.write(f"{a}\t{b}\t{k}\t{len(papers_a)}\t{len(papers_b)}\t{matrix[index[a], index[b]]:.10g}"
                    f"\t{ablated:.10g}\t{empty}\t{np.mean(random_scores):.10g}\t{np.std(random_scores):.10g}"
                    f"\t{capped}\n")


if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    req = aparser.add_argument_group("required arguments")
    req.add_argument("--mesh", help="MeSH descriptor file used for the similarity (e.g. d2025.bin)", required=True)
    req.add_argument("--omim-to-pubmed", help="omim2pubmed file", required=True)
    req.add_argument("--pubmed-to-mesh", help="pubmed2mesh file", required=True)
    req.add_argument("--omim-to-mesh", help="omim2mesh file the similarity was computed from", required=True)
    req.add_argument("--matrix", help="dense similarity matrix (.npy)", required=True)
    req.add_argument("--matrix-ids", help="OMIM ids of the matrix rows, one per line", required=True)
    req.add_argument("--outfile", "--o", help="path to the output file", required=True)
    aparser.add_argument("--categories", default="TWO", choices=CATEGORIES.keys(),
                         help="MeSH categories used for the similarity (default TWO)")
    aparser.add_argument("--random-repeats", type=int, default=20, help="repeats of the random control")
    aparser.add_argument("--validation-sample", type=int, default=20000,
                         help="random non-sharing pairs added to the validation")
    aparser.add_argument("--tolerance", type=float, default=1e-6,
                         help="largest accepted difference between recomputed and matrix scores")
    aparser.add_argument("--seed", type=int, default=0)
    run(aparser.parse_args())
