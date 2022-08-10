"""Provides a utility to extract unique identifiers from the file produced by OMIM_query.py

The output of that script is a mapping with the folowing format:

    OMIM ID\tPubMed ID\tPubMed ID\t...

This script simply collects all pubmeds and writes a file with unique PubMed IDs
"""

import rich.progress

def run(mapping_file, outfile):
    pubmeds = set()
    with rich.progress.open(mapping_file, "r", description="Reading mapping file...") as f:
        for line in f:
            pubmeds |= set(line.strip().split("\t")[1:])
    with open(outfile, "w") as o:
        for pubmed in rich.progress.track(pubmeds, description="Writing unique PubMed IDs..."):
            o.write(f"{pubmed}\n")

if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    required_arguments = aparser.add_argument_group("required arguments")
    required_arguments.add_argument("--omim-to-pubmed", "--in",
                                    help="OMIM2PubMed mapping file",
                                    required=True)
    required_arguments.add_argument("--outfile", "--o",
                                    help="path to the output file", 
                                    required=True)
    args = aparser.parse_args()
    run(args.omim_to_pubmed, args.outfile)

