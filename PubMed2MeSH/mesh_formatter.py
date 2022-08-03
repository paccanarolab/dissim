"""Provide an utility to format a MeSH ASCII file into a simple 2 column file with format

Name\tMeSH ID
"""

import rich.progress

def run(mesh_file, outfile):
    inrecord = False
    outfile = open(outfile, "w")
    title = ""
    identifier = ""
    with rich.progress.open(mesh_file, "r", 
                            encoding="utf8",
                            description="Processing MeSH file...") as f:
        for line in f:
            if "*NEWRECORD" in line:
                inrecord = True
                continue
            if inrecord:
                parts = line.strip().split("=")
                if "MH" == parts[0].strip():
                    title = parts[1].strip()
                if "UI" == parts[0].strip():
                    identifier = parts[1].strip()
                if title and identifier:
                    outfile.write(f"{title}\t{identifier}\n")
                    title = ""
                    identifier = ""
    outfile.close()

if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()

    req = aparser.add_argument_group("required arguments")
    req.add_argument("--mesh-file", "--in",
                     help="MeSH file in ASCII format",
                     required=True)
    req.add_argument("--outfile", "--o", 
                     help="Path to the output file",
                     required=True)
    args = aparser.parse_args()
    run(args.mesh_file, args.outfile)
