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

from collections import defaultdict
from thesaurus import *


class Annotation(object):
    """
    An annotation is a mapping between a set of descriptors (strings)
    and a set of objects (strings). One object can be annotated to several
    descriptors in the thesaurus, and conversely, one descriptor
    in the thesaurus can be annotated to many objects


    For example, the object could be a MIM number.
    """
    def __init__(self):
        self.objects = defaultdict(set)
        self.descriptors = defaultdict(set)
        self.direct_annotations = defaultdict(set)
        self.valid = set()

    def set_valid_descriptors(self, valid):
        """
        Presets a set of valid descriptors. No annotations
        for nodes outside this set will be considered
        """
        self.valid = valid

    def annotate(self, id_object, id_descriptor, node):
        #check if the node is valid
        if (not self.valid) or (id_descriptor in self.valid):
            #add the direct annotation
            self.direct_annotations[id_object].add(id_descriptor)
            for node_anc in node.get_ancestors():
                anc_descriptor = node_anc.get_identifier()
                #check if the ancestor is valid
                if (not self.valid) or (anc_descriptor in self.valid):
                    self.descriptors[anc_descriptor].add(id_object)
                    self.objects[id_object].add(anc_descriptor)

    def get_direct_annotations(self,obj):
        return self.direct_annotations[obj]

    def get_objects(self):
        return self.objects.keys()

    def get_descriptors_per_object(self, obj):
        return self.objects[obj]

    def get_descriptors(self):
        return self.descriptors.keys()

    def num_annot_per_descriptor(self, id_descriptor):
        return len(self.descriptors[id_descriptor])

    def get_objects_per_descriptor(self,descriptor):
        """
        Returns the diseaes that this node annotates
        """
        return self.descriptors[descriptor]

class AnnotationParser(object):
    """
    Creates an annotation from a file and a Thesaurus,
    where each row is an object (first column) and the
    rest of the fields are descriptors. Annotations are
    up-propagated following the true path rule (when a
    node is annotated, also are its ancestors)

    the file is the file formed by:
        OMIM\tMESH_1\tMESH_2..MESH_K.
    
    """
    def __init__(self, thesaurus, datafile):
        self.thesaurus = thesaurus
        self.__data = defaultdict(list)
        #read the annotation file.
        self.__readAnnotations(datafile)

    def __readAnnotations(self, datafile):
        with open(datafile) as f:
            for line in f:
                fields = line.strip().split()
                self.__data[fields[0]] = fields[1:]


    def get_annotations(self, chosen_categories=[]):
        annot = Annotation()
        descriptors_in_thesaurus = set()
        #if something is specified, fetch only these nodes.
        if chosen_categories:
            for cat in chosen_categories:
                descriptors_in_thesaurus.update(self.thesaurus.get_nodes_by_category(cat))
            annot.set_valid_descriptors(descriptors_in_thesaurus)
        else: #nothing is speciefied, bring everything
            descriptors_in_thesaurus = set(self.thesaurus.get_node_ids())

        #read the annotation file
        for obj in self.__data:
            id_descriptors = set(self.__data[obj]) & descriptors_in_thesaurus
            for id_descriptor in id_descriptors: 
                node = self.thesaurus.get_node(id_descriptor)
                annot.annotate(obj,id_descriptor,node)
        return annot
