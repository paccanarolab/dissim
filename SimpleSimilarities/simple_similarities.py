
"""
    Calculates simple semantic similarity measures
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
import itertools as it
import sys
import progressbar

#descriptor-tree positions file
def readTreePositions(filename):
    node_coord = defaultdict(list)
    with open(filename,'r') as f:
        for line in f:
            sl = line.strip().split()
            node_coord[sl[0]] = set([i.strip().split('.')[0][0] for i in sl[1::]])
    return node_coord

#mim2mesh
def read_mapping(filename):
    d = defaultdict(set)
    with open (filename, "r") as f:
        for line in f:
            sl = line.strip().split()
            d[sl[0]] = set(sl[1:])
    return d

def filterAnnotation(mapping,desired_trees,tree_positions):
    filtered_mapping = defaultdict(list)
    for val in mapping:
        for desc in mapping[val]:
            if set(desired_trees) & set(tree_positions[desc]):
                filtered_mapping[val].append(desc)

    return filtered_mapping

def jaccard(set1, set2):
    try:
        return float(len(set1 & set2)) / float(len(set1 | set2))
    except:
        return 0
        
def dice(set1, set2):
    try:
        return float(2 * len(set1 & set2)) / float(len(set1) + len(set2))
    except:
        return 0

def overlap(set1, set2):
    try:
        return float(len(set1 & set2)) / float(min(len(set1), len(set2)))
    except:
        return 0

def num_common(set1, set2):
    return float(len(set1 & set2))

def compute_combined(values, output_folder, filename_modifier, method_key,tree_positions, filter_list):
    options = {0 : jaccard, 1 : dice, 2: overlap, 3 : num_common}
    names = {0 : "jaccard", 1 : "dice", 2: "overlap", 3 : "num_common"}
    print(names[method_key])
    filtered_values =  filterAnnotation(values,filter_list, tree_positions)
    name = output_folder + filename_modifier + '_' + names[method_key]
    #merge them
    print('Calculating ' + name)
    with open(name, "w") as out:
        for d1, d2 in it.combinations(sorted(values.keys()), 2):
            s1 = set(filtered_values[d1])
            s2 = set(filtered_values[d2])
            similarity = options[method_key](s1,s2)
            if similarity > 0.0:
                out.write("%s\t%s\t%f\n" % (d1, d2, similarity))

def compute_per_ontology(values,output_folder, key, tree_positions):
    options = {0 : jaccard, 1 : dice, 2: overlap, 3 : num_common}
    names = {0 : "jaccard", 1 : "dice", 2: "overlap", 3 : "num_common"}

    available_ontologies = dict([
        ('A','Anatomy'), ('B','Organisms'),
        ('C','Diseases'), ('D','Chemicals and Drugs'),
        ('E','Analytical, Diagnostic and Therapeutic Techniques and Equipment'),
        ('F','Psychiatry and Psychology'),
        ('G','Phenomena and Processes'),
        ('H','Disciplines and Occupations'),
        ('I','Anthropology, Education, Sociology and Social Phenomena'),
        ('J','Technology, Industry, Agriculture'),
        ('K','Humanities'),
        ('L','Information Science'), ('M','Named Groups'),
        ('N','Health Care'), ('Z','Geographicals')
    ])
    for current_ontology in available_ontologies:
        name = output_folder + current_ontology + '_' + names[key]
        filtered_values =  filterAnnotation(values, current_ontology, tree_positions)
        #bar = progressbar.ProgressBar(maxval = (len(filtered_values.keys())*(len(filtered_values.keys()) - 1))/2, widgets=[progressbar.Bar('*', '[', ']'), ' ', progressbar.SimpleProgress()]).start()
        #merge them
        barCounter = 0
        zero = 0
        total = 0
        print('Calculating ' + name)
        with open(name, "w") as out:
            for d1, d2 in it.combinations(sorted(filtered_values.keys()), 2):
                s1 = set(filtered_values[d1])
                s2 = set(filtered_values[d2])
                similarity = options[key](s1,s2)
                if similarity > 0.0:
                    out.write("%s\t%s\t%f\n" % (d1, d2, similarity))
                else:
                    zero += 1
                total += 1
             #   bar.update(barCounter)
             #   barCounter+=1

        print(current_ontology + ' zero: ' + str(zero) + ' total: '+  str(total))

help_string = """
        ---------------------------------------------------------------------------------------------------------------
        Computes simple similarities based on shared elements (MeSH terms/Publications) for the diseases
        Usage: 
        =======
        sys.argv[0] mapping_file descriptor output_folder
        \t* mapping_file: Either and OMIM to PubMed mapping or an OMIM to MeSH mapping. Produced by MIM2MESH.py or Pubmed_query.py
        \t* descriptor_tree_position: Mapping between MeSH descriptors and their tree coordinates. 
            The file format: MeSHDescriptor\\tCoord1\\tCoord2...\\tCoordN
        \t* output_folder: Folder where the results will be placed.
        ---------------------------------------------------------------------------------------------------------------
        """


if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print(help_string)
        sys.exit(-1)

    output_folder = sys.argv[3]

    tree_positions = readTreePositions(sys.argv[2])
    values = read_mapping(sys.argv[1])
    
    print('Per ontology ',)
    for i in range(0,4):
        compute_per_ontology(values, sys.argv[3], i,tree_positions)

    print('Combined ',)
    allontologies = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Z']
    fiveontologies =  ['A','C','D','E','G']
    for i in range(0,4):
        compute_combined(values, output_folder, '5categories', i ,tree_positions, fiveontologies)
        compute_combined(values, output_folder, 'combinedCategories', i ,tree_positions, allontologies)

