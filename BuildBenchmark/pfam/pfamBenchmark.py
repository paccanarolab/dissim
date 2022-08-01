from collections import defaultdict
from itertools import combinations
import sys

"""
    Produce the pfam benchmark
    Copyright (C) 2015 Horacio Caniza

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

__author__  = "Horacio Caniza"
__email__   = "h.j.canizavierci@cs.rhul.ac.uk"
__copyright__ = "Copyright (C) 2015 Horacio Caniza"
__license__ = "GPL"
__version__ = "3"


def readPfamScan(filename):
    values = defaultdict(set)
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == "#" or line[0] == "":
                continue
            sl = line.strip().split()
            values[sl[0].split('|')[1]].add(sl[6].strip().upper())
    return values

def readMimToProt(filename):
    values = defaultdict(list)
    with open(filename, 'r') as f:
        for line in f:
            sl = line.strip().split()
            values[sl[0]].append(sl[1])
    return values

def readExclusionList(filename):
    values = set()
    with open(filename, 'r') as f:
        for line in f:
            values.add(line.strip().upper())
    return values

def produceBenchmark(pfamscan, mimtoprot, exclusion, outfile_name):
    valid_pairs = list()
    with open(outfile_name,'w') as f:
        for (omim1, omim2) in combinations(sorted(mimtoprot.keys()), 2):
            #get the proteins
            proteins_omim1 = mimtoprot[omim1]
            #get the families
            families_omim1 = set()
            for protein in proteins_omim1:
                families_omim1.update(pfamscan[protein])
            #get the proteins for the second omim
            proteins_omim2 = mimtoprot[omim2]
            #get the families
            families_omim2 = set()
            for protein in proteins_omim2:
                families_omim2.update(pfamscan[protein])
            #get the common families

            common_families = set(families_omim1) & set(families_omim2)
            #if the exclusion is empty, then we can add a one.
            if common_families - exclusion:
                f.write(str(omim1) + '\t' + str(omim2) + "\t1\n")


help_str = """
Produces a binary disease similarity graph
Use: python pfamBenchmark.py pfamscan omim_to_uniprot outfile_name [exclusion]'

It takes two files as an input:
     1.- The output of pfam_scan on a set of protein
     sequences 

     2.- a table relating omim diseases and
     protein ids (one entry per pair)
         omim_i    protein_j
         omim_i    protein_k
         ...
  An additional argument is a single column file
  that lists the pfams to be excluded

The third parameter will be an output file 
with the format 
         omim_i    omim_j    1.0

if two diseases have proteins associated shared 
at least one pfam family
"""

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(help_str)
        exit()

    pfamscan_filename = sys.argv[1]
    mimtoprot_filename = sys.argv[2]
    outfile_name = sys.argv[3]
    
    values_exclusion = set()
    if len(sys.argv) == 5:
        exclusionList_filename = sys.argv[4]
        print('Reading exclusion list..')
        values_exclusion = readExclusionList(exclusionList_filename)
    
    print('Reading pfamscan file..')
    values_pfamscan = readPfamScan(pfamscan_filename)
    print('Reading mimtoprot file..')
    values_mimtoprot = readMimToProt(mimtoprot_filename)
    #just to see if the fourt parameter was provided
    print('Producing benchmark..')
    produceBenchmark(values_pfamscan, values_mimtoprot, values_exclusion, outfile_name)
