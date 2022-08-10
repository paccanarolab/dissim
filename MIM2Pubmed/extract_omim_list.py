"""Provides a utility to extract omim identifiers from the mimTitles.txt file provided by OMIM

This is mainly provided for people working in Windows environments, 
if you have access to cut and grep, simply run:

cut -f2 mimTitles.txt | grep -v '#' > omim_ids.txt 

and you will have the same output as this script
"""

import rich.progress

def run(mim_titles, outfile):
    omim_ids = set()
    with rich.progress.open(mim_titles, "r", description="Extracting OMIM IDs...") as f:
        for line in f:
            if line[0] == "#":
                continue
            omim_ids.add(line.strip().split("\t")[1])
    with open(outfile, "w") as o:
        for omim_id in rich.progress.track(omim_ids, description="Writing output..."):
            o.write(f"{omim_id}\n")

if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    req = aparser.add_argument_group("required arguments")
    req.add_argument("--mim-titles", "--in", 
                     help="path to mimTitles.txt by OMIM",
                     required=True)
    req.add_argument("--outfile", "--o",
                     help="path to the output file", 
                     required=True)
    args = aparser.parse_args()
    run(args.mim_titles, args.outfile)

