"""
    Computes semantic similarity in the MeSH ontologies.
    Copyright (C) 2015 Alfonso E. Romero, Horacio Caniza

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

__author__  = "Alfonso E. Romero, Horacio Caniza"
__email__   = "h.j.canizavierci@cs.rhul.ac.uk"
__copyright__ = "Copyright (C) 2015 Horacio Caniza"
__license__ = "GPL"
__version__ = "3"

from thesaurus import *
from annotation import *
import numpy as np
from collections import OrderedDict
from rich.progress import track

class SemanticSimilarity(object):

    def __init__(self, thesaurus, annotation):
        self.annotation = annotation
        self.thesaurus = thesaurus
        self.objects = list(self.annotation.get_objects())
        self.descriptors = list(self.annotation.get_descriptors())
        self.num_objects = len(self.objects)
        self.num_descriptors = len(self.descriptors)
        self.lowestCommonAncestor = defaultdict(list)
        self.selectedPair = defaultdict(list)

        self.object_indexes = OrderedDict()
        self.descriptors_indexes = OrderedDict()
        for i,desc in enumerate(self.descriptors):
            self.descriptors_indexes[desc] = i
        for i,obj in enumerate(self.objects):
            self.object_indexes[obj] = i


    def semantic_similarity(self, id1, id2):
        raise NotImplementedError("SemanticSimilarity is an abstract class. Check similarity_measures.py for the correct class implementation")

    def get_descriptor_indexes(self):
        return self.descriptors_indexes

    def get_object_indexes(self):
        return self.object_indexes

    def get_perDescriptor(self):
        try:
            return self.perDescriptor
        except:
            raise

    def get_perObject(self):
        try:
            return self.perObject
        except:
            raise

    def get_selectedPair(self):
        return self.selectedPair

    def get_lowestCommonAncestor(self):
        return self.lowestCommonAncestor


    def common_ancestors(self, id1, id2):
        return [node for node in self.thesaurus.common_ancestors(id1, id2) if node.get_identifier() in self.descriptors]
            
    def __get_descriptors_ids_per_object(self):

        descriptors_id_per_object = defaultdict(set)
        for i in range(self.num_objects):
            descriptors_id_per_object[i] = [self.descriptors_indexes[j] for j in self.annotation.get_direct_annotations(self.objects[i])]

        return descriptors_id_per_object



    def compute_semantic_similarity_per_object_diseasewise(self):
        self.perObject = np.zeros((self.num_objects, self.num_objects))
        for dis1 in track(range(self.num_objects), description="Computing semantic similarity..."):
            for dis2 in range(dis1, self.num_objects):
                similarity = self.semantic_similarity(self.objects[dis1], self.objects[dis2])
                self.perObject[dis1, dis2] = similarity
                self.perObject[dis2, dis1] = similarity
        

    def compute_semantic_similarity_per_object_termwise(self):

        self.perObject = np.zeros((self.num_objects, self.num_objects))
        descriptors_id_per_object = self.__get_descriptors_ids_per_object()
        with open('./localStore/sim_distribution','w') as f:
            f.write('#disease A\tdisease B\tmax_similarity\tmin_similarity\tmean_similarity\tmedian_similarity\tstandard_deviation\n')
            for i in track(range(self.num_objects), description="Computing term-wise similarity..."):
                rows = self.perDescriptor[ descriptors_id_per_object[i] ]
                for j in range(i, self.num_objects):
                    similarity = self.selectionStrategy(rows[:,descriptors_id_per_object[j]])

                    submat = rows[:,descriptors_id_per_object[j]]
                    mean = np.mean(submat)
                    median = np.median(submat)
                    std_dev = np.std(submat)
                    f.write(self.objects[i] + '\t' + self.objects[j] + '\t' + str(np.max(submat)) + '\t' + str(np.min(submat)) + '\t' + str(mean) + '\t' + str(median) + '\t' + str(std_dev)+'\n')
                    self.perObject[i, j] = similarity
                    self.perObject[j, i] = similarity


    def compute_semantic_similarity_per_descriptor(self):

        self.perDescriptor = np.zeros((self.num_descriptors, self.num_descriptors))
        max_dist = -1000
        for i in track(range(self.num_descriptors), description="Computing semantic similarity per descriptor..."):
            desc1 = self.descriptors[i]
            for j in range(i, self.num_descriptors):
                desc2 = self.descriptors[j]
                (selectedAncestor, similarity) = self.semantic_similarity(desc1, desc2)
                #assign the similarities.
                self.perDescriptor[i, j] = similarity
                self.perDescriptor[j, i] = similarity
                #store the node that was selected. Quicker to store both, as speed is mor of a problem
                self.lowestCommonAncestor[desc1].append((desc2,selectedAncestor,similarity))

        #we need to call a normalisation function because of jiang. Each measure 
        #implemets their own normalisation if needed.
        self.normalise(self.perDescriptor)


