"""Provides a utility to extract omim identifiers from the files provided by OMIM

Two input files are supported, and the format is detected from the file header:

* mim2gene.txt: freely available, no registration needed
  (https://omim.org/static/omim/data/mim2gene.txt). The MIM number is the first column.
* mimTitles.txt: requires OMIM download registration. The MIM number is the second column.

Moved/removed entries (type "moved/removed" in mim2gene.txt, prefix "Caret" in mimTitles.txt)
are skipped by default, as they do not reference any publication. Use --keep-removed to keep them.

This is mainly provided for people working in Windows environments,
if you have access to grep and awk, simply run:

grep -v '^#' mim2gene.txt | awk -F'\\t' '$2 != "moved/removed" {print $1}' > omim_ids.txt

and you will have the same output as this script
"""

import rich.progress

MIM2GENE = "mim2gene"
MIM_TITLES = "mimTitles"


def detect_format(header_lines):
    for line in header_lines:
        if line.startswith("# MIM Number"):
            return MIM2GENE
        if line.startswith("# Prefix"):
            return MIM_TITLES
    raise ValueError("Could not detect the file format, expected mim2gene.txt or mimTitles.txt")


def run(omim_file, outfile, keep_removed=False):
    header_lines = []
    rows = []
    with rich.progress.open(omim_file, "r", description="Extracting OMIM IDs...") as f:
        for line in f:
            if line[0] == "#":
                header_lines.append(line)
                continue
            rows.append(line.rstrip("\n").split("\t"))
    file_format = detect_format(header_lines)

    omim_ids = set()
    for fields in rows:
        if file_format == MIM2GENE:
            mim, removed = fields[0], fields[1] == "moved/removed"
        else:
            mim, removed = fields[1], fields[0] == "Caret"
        if removed and not keep_removed:
            continue
        omim_ids.add(mim.strip())
    with open(outfile, "w") as o:
        for omim_id in rich.progress.track(sorted(omim_ids), description="Writing output..."):
            o.write(f"{omim_id}\n")

if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    req = aparser.add_argument_group("required arguments")
    req.add_argument("--omim-file", "--mim-titles", "--in",
                     help="path to mim2gene.txt or mimTitles.txt by OMIM",
                     required=True)
    req.add_argument("--outfile", "--o",
                     help="path to the output file",
                     required=True)
    aparser.add_argument("--keep-removed", action="store_true",
                         help="keep moved/removed entries (they have no references)")
    args = aparser.parse_args()
    run(args.omim_file, args.outfile, args.keep_removed)
