#!/usr/bin/python

import sys
from collections import defaultdict
import progressbar

def readMappingFile(filename):
    values = defaultdict(list)
    with open(filename, 'r') as mappingFile:
        for line in mappingFile:
            sl = line.strip().split()
            values[sl[0]] = sl[1:]
    return values

def mim2mesh(mim2pubmed, pubmed2mesh, outfile):
    bar = progressbar.ProgressBar(maxval = len(mim2pubmed), widgets=[progressbar.Bar('*', '[', ']'), ' ', progressbar.Percentage()]).start()
    #merge them
    barCounter = 0
    with open(outfile,'w') as f:
        for mimno in mim2pubmed:
            line = [mimno]
            for pubmed in mim2pubmed[mimno]:
                line.extend(pubmed2mesh[pubmed])
            if len(line) > 1:
                f.write('\t'.join(line))
                f.write('\n')

            #update the progressbar
            barCounter = barCounter + 1
            bar.update(barCounter)

    bar.finish()

if __name__ == '__main__':
    
    #check the parameters.
    if len(sys.argv[1:]) < 3:
        print("--------------------------------------------------------------------------------------------------------------")
        print("Maps each OMIM record to the MeSH terms of its associated publication")
        print("Usage:")
        print("======")
        print("python MIM2MeSH.py mim2pubmed pubmed2mesh outputfile")
        print("Where:")
        print("mim2pubmed: Mapping between OMIM records and the PubMed identifiers of the referenced literature ")
        print("\t*mim\\tpubmedid\\t...\\tpubmedid")
        print("pubmed2mesh: Mapping between PubMed records and their MeSH terms")
        print("\t*pubmedid\\tmeshDescriptorUniqueId\\t...\\tmeshDescriptorUniqueId")
        print("outputfile: Desired output file name")
        print("--------------------------------------------------------------------------------------------------------------")
        exit()
    #------------------

    #read the first file.
    mim2pubmed = readMappingFile(sys.argv[1])
    #read the second file.
    pubmed2mesh = readMappingFile(sys.argv[2])
    #produce the mapping
    mim2mesh(mim2pubmed, pubmed2mesh, sys.argv[3])
