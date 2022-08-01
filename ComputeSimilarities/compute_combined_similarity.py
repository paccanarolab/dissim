"""
    Computes semantic similarity in the MeSH ontologies.
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


help_string = """
    --------------------------------------------------------------------------------------------------------------
    Usage:
    =====
    python compute_combined_similarity.py descriptors_file annotation_file chosen_measure ism category_subset [all/two/five] Optional:filename_modifier

    *Small format guide:

        descriptors_file:
            The ASCII version of MeSH.
            Format:
                See d2014.bin file provided by MeSH.


        annotation_file:
            File produced by MIM2MESH.py
            Format:
                OMIM\\tMeSH_descriptorUID...\\tMeSH_descriptorUID

        chosen_measure:
            Determines the measure that will be calculated
                Values:
                    Resnik
                    Lin
                    Jiang
                    SimUI
                    SimGIC

        ism:
            Determines whether the ISM will be computed for the selected measure. See:
            Improving GO semantic similarity measures by exploring the ontology beneath the terms and modelling uncertainty. Yang et al.
            For the results provided in the paper, this was set to no.
                Values:
                    no (default)
                    yes

        category_subset:
            The subset of MeSH ontologies that will be combined to produce the disease similarity data.  
                Values:
                    All: All MeSH ontologies
                    Two: A and C ontologies
                    Five: A, C, D, E and G ontologies

        filename_modifier:
                Appends a string modifier to the filename, for flexilibity mainly. If nothing is set, that's ok.

    --------------------------------------------------------------------------------------------------------------

"""


import sys
import os
from mesh_parser import MeSHParser
from thesaurus import *
from annotation import *
from semsim import *
from similarity_measures import *
from writeFiles import *
import progressbar

#CATEGORIES
categories = dict([
    ('TWO',['A','C']),
    ('FIVE',['A','C','D','E','G']),
    ('ALL',['A','B','C','D','E','F','G','H','I','J','K','L','M','N','V','Z'])
    ])

methods_termwise = {0: (Resnik,'MAX'), 1:(Lin,'MAX'), 2: (Jiang,'MAX'), 3: (Schlicker, 'MAX')}
names_termwise = {"RESNIK":0, "LIN":1, "JIANG":2, "SCHLICKER":3}

methods_diseasewise = {0: SimUI, 1:SimGIC}
names_diseasewise = {"SIMUI":0, "SIMGIC":1}



if len(sys.argv) < 5:
    print help_string
    sys.exit(-1)

#0. prepare utilities

descriptors_file = sys.argv[1]
annotation_file = sys.argv[2]
chosen_measure = sys.argv[3].upper()
compute_ism = sys.argv[4].upper()
categories_subset = sys.argv[5].upper()

file_per_disease = 'combined_similarity-' + categories_subset + "_"+ chosen_measure
if len(sys.argv) == 7:
    file_per_disease = file_per_disease + sys.argv[6]

# 1.- we load up the thesaurus with the categories we want.
parser = MeSHParser(descriptors_file, categories[categories_subset])
thesaurus = parser.get_thesaurus(True)

print "\t- Obtaining annotation"
annotation_parser = AnnotationParser(thesaurus, annotation_file)
annotation = annotation_parser.get_annotations()

print "\t- Computing " + chosen_measure
if chosen_measure in names_termwise:
    method_pair = methods_termwise[names_termwise[chosen_measure]]
    sem_sim = method_pair[0](thesaurus, annotation,method_pair[1])
    #---------
    #per decriptor
    print '\t\t- Computing per descriptor..'
    sem_sim.compute_semantic_similarity_per_descriptor()
    #per object
    print '\t\t- Computing per object..'
    sem_sim.compute_semantic_similarity_per_object_termwise()
    #--------------
    ##per descriptor
    cache_file = './Cache/combined_' + chosen_measure + '_' + str(categories_subset) +"_per_descriptor.txt"
    if os.path.isfile(cache_file):
        sem_sim.perDescriptor = np.loadtxt(cache_file,delimiter='\t')
    else:
        print '\t\t- Computing per descriptor..'
        sem_sim.compute_semantic_similarity_per_descriptor()
        print '\t\t- Writing per descriptor'
        np.savetxt(cache_file ,sem_sim.get_perDescriptor(),delimiter='\t', newline='\n')
    #per object
    cache_file = './Cache/'+ chosen_measure + "_combined_" +str(categories_subset) +"_per_disease.txt"
    if os.path.isfile(cache_file):
        sem_sim.perObject = np.loadtxt(cache_file,delimiter='\t')
    else:
        print '\t\t- Computing per object..'
        sem_sim.compute_semantic_similarity_per_object_termwise()
        print '\t\t- Saving per disease'
        np.savetxt(cache_file ,sem_sim.get_perObject(),delimiter='\t', newline='\n')
    ##--------------
    print '\t\t-Get LCA..'
    lowest_common_ancestor = sem_sim.get_lowestCommonAncestor()
elif chosen_measure in names_diseasewise:
    sem_sim = methods_diseasewise[names_diseasewise[chosen_measure]](thesaurus,annotation)
    print '\t- Calculating per disease..'
    sem_sim.compute_semantic_similarity_per_object_diseasewise()

per_disease = sem_sim.get_perObject()

print "\t -Writing file.."

writeTriplet(file_per_disease, per_disease, sem_sim)

print "\t -Write LCA..."
writeSelectedDescriptor(file_per_disease +"-LCA", lowest_common_ancestor)

if compute_ism == "YES":
    print "\t -Computing ISM for " + chosen_measure
    ism = ISM(thesaurus, annotation, per_disease)
    ism.ism()
    per_disease= ism.ISM
    file_per_disease = 'ISM-combined_similarity-' + chosen_measure
    print "\t -Writing file for ISM_"+chosen_measure
    writeTriplet(file_per_disease, per_disease, sem_sim)
