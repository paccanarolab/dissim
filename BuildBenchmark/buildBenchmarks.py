
"""
    Builds the benchmarks
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


import sys
import tempfile
from scipy.stats.stats import pearsonr
import itertools

from os import listdir
from os.path import isfile, join

def readDiseasesWithProteins(filename):
    values = set()
    with open(filename, 'r') as f:
        for line in f:
            values.add(int(line.strip()))
    return values

def get_ids_from_file(filename):
    """
    Reads a matrix file and gets in a set
    all the elements in the first two columns
    (i.e. OMIM diseases)
    """
    ids_in_file = set()
    with open(filename, 'r') as f:
        for line in f:
            (omim1, omim2, _) = line.split()
            ids_in_file.add(int(omim1))
            ids_in_file.add(int(omim2))
    return ids_in_file

def extract_subset(filename, ids):
    matrix = {}
    with open(filename, 'r') as f:
        for line in f:
            (omim1, omim2, simi) = line.strip().split()
            omim1 = int(omim1)
            omim2 = int(omim2)
            if omim1 in ids and omim2 in ids and omim1 != omim2:
                mini, maxi = sorted([omim1, omim2])
                matrix[(mini,maxi)] = float(simi) 
    return matrix

 
def print_intersection_matrix(file1, file2, ids, output):
    """
    Take the matrix file `file1` and the matrix file
    `file2`, the intersection set of OMIM ids `ids`
    and the name of the output file, and prints an
    output file
        OMIM_ID_1   OMIM_ID_2   SIM_1   SIM2 
        ...
    where for every pair of OMIM disease ids in the
    intersection (columns 1 and 2) prints to the output
    file the similarities in the first file and in 
    the second file in columns 3 and 4, respectively.

    For our proposed measure, we need to satisfy three criteria:

    """
    matrix_ground_truth = extract_subset(file1, ids)
    matrix_disim = extract_subset(file2, ids)

    with open(output, 'w') as out:
        list_ids = sorted(list(ids))
        for (id1, id2) in itertools.combinations(list_ids, 2):
            simi_disim = matrix_disim[(id1, id2)] if (id1, id2) in matrix_disim else 0.0
            simi_ground_truth = matrix_ground_truth[(id1, id2)] if (id1, id2) in matrix_ground_truth else 0.0
            out.write("%i\t%i\t%f\t%f\n" % (id1, id2, simi_disim, simi_ground_truth))
                    

help_string = """ 
        Use: python buildBenchmarks.py  ground_truth_file similarity_file filtered_mimtoprot destination_dir filename_modifier

            ground_truth_file: Molecular similarity file. Format:
                OMIM_ID_1\\tOMIM_ID_2\\t1/0

            similarity_file: Disease similarity file produced by compute_combined_similarities.py, compute_matrices.py, simple_similarities.py. Format:
                OMIM_ID_1\\tOMIM_ID_2\\tSIM

            filtered_mimtoprot: File mapping OMIM diseases to UniProt Indentifiers. Format:
                OMIM1\\tUNIPROT1
                OMIM1\\tUNIPROT2
                OMIM2\\tUNIPROT2
            
            destination_dir: Location of the result files.
            filename_modifier: String appended to the end of the file.
    """


if __name__ == "__main__":
    if (len(sys.argv) != 6):
        print(help_string)
        sys.exit(-1)

    #get the ids in the similarity file.
    ids_in_2 = get_ids_from_file(sys.argv[2])
    print("Number of OMIM ids in ", sys.argv[2], ": ", len(ids_in_2))
    #get the diseaes with proteins.
    diseases_with_proteins = readDiseasesWithProteins(sys.argv[3])
    #compute the intersection.
    testeable_omims = ids_in_2 & diseases_with_proteins
    #print "Common OMIM ids " + str(len(intersection))
    print_intersection_matrix(sys.argv[1], sys.argv[2], testeable_omims, sys.argv[4] + sys.argv[2].split('/')[-1] + '_' + sys.argv[5])
