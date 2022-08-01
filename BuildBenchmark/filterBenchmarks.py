
"""
    Filters the benchmarks.
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

from collections import defaultdict
import sys


def readUniprotPfamMapping(filename, invalid_pfams):
    values = set()
    with open(filename, 'r') as f:
        for line in f:
            #skip comments
            if line[0] == "#" or line[0] == "":
                continue
            #get the different fields.
            sl = line.strip().split()
            #check if the pfam is not invalid.
            if sl[6].strip().upper() not in invalid_pfams:
                values.add(sl[0].split('|')[1])
    return values


def readInvalidPfam(filename):
    values = set()
    with open(filename, 'r') as f:
        for line in f:
            values.add(line.strip())
    return values


def readOmimToEntity(filename):
    values = defaultdict(set)
    with open(filename, 'r') as f:
        for line in f:
            (omim, entity) = line.strip().split()
            values[omim].add(entity)
    return values


def writeLine(fp, omim1, omim2, disim, molsim):
    fp.write("%s\t%s\t%s\t%s\n" % ( omim1, omim2, disim, molsim))

def filterFile(filename, omim_uniprot, valid_uniprot, outfile_name):
    with open(filename,'r') as infile, open(outfile_name, 'w') as outfile:
        for line in infile:
            #get the values in the file to filter.
            (omim1, omim2, disim, molsim) = line.strip().split()
            #filter proteins
            #if not valid_uniprot, this is nothign to do with pfam, so we just check
            #if the pair shares a protein.
            if not valid_uniprot:
                if not (omim_uniprot[omim1] & omim_uniprot[omim2]):
                    writeLine(outfile, omim1,omim2,disim,molsim)
            #check if omim_pfam is true, if it is, then filter.
            #should there be anything in valid pfam, we need to check if they do not share a protein
            #and both are valid uniprots, i.e. do not have a forbidden pfam.
            elif valid_uniprot:
                #since omim_pfam contains the valid pfams for each disease, whenever a disease appears, it is valid.
                if (omim_uniprot[omim1] & valid_uniprot) and (omim_uniprot[omim2] & valid_uniprot) and not (omim_uniprot[omim1] & omim_uniprot[omim2]):
                    writeLine(omim1,omim2,disim,molsim)
            #we are done!

help_str = """
Use: python filterBenchmarks.py file_to_filter omim_uniprot [uniprot_pfam_table invalid_pfams]

    1. file_to_filter: File produced by buildBenchmarks.py

    2. omim_uniprot: Double column file. This file contains, in the first column and OMIM number and in the second column 
    the uniprot id of the proteins associated to the diseaes. Should a disease have more than one protein, it will be repeated
    in as many lines as proteins it has.

    3. uniprot_pfam_table: result of pfamscan.pl. This file contains a mapping of uniprot ids and their pfam families, 
    domains, motifs and repeats. 

    4. invalid_pfams: single column file. Pfams to be excluded because they show up in mesh.

"""

if __name__ == "__main__":
    #if we have one optional parameter, but not the second one.
    if len(sys.argv) < 3 or (len(sys.argv) > 3 and len(sys.argv) < 5):
        print(help_str)
        exit()
    #set the input parametesr.
    file_to_filter = sys.argv[1]
    omim_uniprot_file = sys.argv[2]
    uniprot_pfam_file = ""
    invalid_pfam_file = ""

    if len(sys.argv) == 5:
        uniprot_pfam_file = sys.argv[3]
        invalid_pfam_file = sys.argv[4]

    #read omim to uniprot mapping
    print('Reading omim_to_uniprot..')
    omim_uniprot_mapping = readOmimToEntity(omim_uniprot_file)
    #check if omim to pfam is required.
    valid_uniprots_pfam = set()
    invalid_pfam = []
    if uniprot_pfam_file:
        print('Reading invalid pfams..')
        invalid_pfam = readInvalidPfam(invalid_pfam_file)
        print('Reading uniprot to pfam mapping..')
        valid_uniprots = readUniprotPfamMapping(uniprot_pfam_file, invalid_pfam)

    #read the file to filter, considering the filters and produceoutput.
    print('Filtering file..')
    filterFile(file_to_filter, omim_uniprot_mapping, valid_uniprots_pfam, file_to_filter+'-filter')
