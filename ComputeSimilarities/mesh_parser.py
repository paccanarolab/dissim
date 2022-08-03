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

import sys
from collections import defaultdict
from thesaurus import MeSHThesaurus
from thesaurus import MeSHThesaurusNode
import re
from itertools import product
from annotation import *
import numpy as np


""" 
    A simple MeSH thesaurus parser, Requires:
    - the ASCII version of MeSH (which can be downloaded
        from the site). Descriptor file d2013.bin
    The get_thesaurus() method will invoke the parser and
    a MeSH thesaurus object will be returned
"""


class MeSHParser(object):
    def __init__(self, mesh_file,categories):
        #the categories in MeSH
        self.available_categories = dict([
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
            ('N','Health Care'), ('V','Publication Characteristics'), ('Z','Geographicals')
        ])
        self.categories = categories
        self.mesh_file = mesh_file

    # 1.- we create the nodes and add them to the thesaurus
    #it will create nodes at level 0, that is nodes corresponding to :
    #Anantomy,Disciplines, etc. So the position will be: A, D, etc.
    def __create_root_nodes(self, thesaurus):
        for cat in self.categories:
            node = MeSHThesaurusNode(self.available_categories[cat], cat, [cat])
            node.is_dummy = True
            thesaurus.add_category(cat, self.available_categories[cat])
            thesaurus.add_node(node)


    def __add_generic_node(self, thesaurus):
            generic_node = 'GEN'
            node = MeSHThesaurusNode(generic_node, generic_node, [generic_node])
            node.is_dummy = True
            thesaurus.add_category(generic_node, generic_node)
            thesaurus.add_node(node)
            #link each category to the new generic root node.
            for category_id in thesaurus.get_category_ids():
                category_node = thesaurus.get_node(category_id)
                category_node.add_parent(node)
                node.add_child(category_node)


    def __parse_MeSH_file(self, thesaurus):
        r = re.compile("\\s+")
        # 1.- we parse the MeSH descriptor file
        with open(self.mesh_file, encoding="utf8") as f:
            tree_positions = []
            synonyms = []
            name = ""
            identifier = ""
            is_a_tree = False
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    #we need to check if the current node is supposed to be added
                    #based on the ontologies we have selected.
                    tree_positions = [i for i in tree_positions if i[0] in thesaurus.get_category_ids()]
                    if tree_positions:
                        node = MeSHThesaurusNode(name, identifier, tree_positions)
                        node.set_synonyms(synonyms)
                        if is_a_tree:
                            node.set_dummy(True)
                            thesaurus.add_tree(identifier, name)
                        thesaurus.add_node(node)

                    #restart the whole process.
                    name = ""
                    identifier = ""
                    is_a_tree = False
                    tree_positions = []
                    synonyms = []
                elif line.startswith("UI = "):
                    identifier = line.split()[2].strip()
                elif line.startswith("MN = "):
                    tree_pos = line.split()[2].strip()
                    is_a_tree = is_a_tree or tree_pos.find('.') == -1
                    tree_positions.append(tree_pos)
                elif line.startswith("MH = "):
                    name = line.split('=')[1].strip()
                elif line.startswith("ENTRY = ") or line.startswith("PRINT ENTRY = "):
                    l = line.split("=", 1)[1].strip()
                    if "|" in l:
                        l = l.split("|")[0]
                        synonyms.append(l)

                
            if len(tree_positions) != 0:
                node = MeSHThesaurusNode(name, identifier, tree_positions)
                if is_a_tree:
                    node.set_dummy(True)
                    thesaurus.add_tree(identifier, name)
                thesaurus.add_node(node)
    
        # 2.- we adjust relationships between trees and categories
        for tree_id in thesaurus.get_trees_ids():
            tree_node = thesaurus.get_node(tree_id)
            #several categories even for the "root" of the original trees.
            for category_id in tree_node.get_categories():
                category_node = thesaurus.get_node(category_id)
                category_node.add_child(tree_node)
                tree_node.add_parent(category_node)

        # 3.- we adjust the rest of relationships
        for node_id in thesaurus.get_node_ids():
            node = thesaurus.get_node(node_id)
            if not node.is_dummy:
                for tree_position in node.get_tree_positions():
                    parent_position = '.'.join(tree_position.split('.')[0:-1]).strip()
                    parent = thesaurus.get_node_by_position(parent_position)
                    parent.add_child(node)
                    node.add_parent(parent)
            for cat_id in node.get_categories():
                thesaurus.nodes_by_category[cat_id].add(node_id)
        ####---

    def get_thesaurus(self, with_root = False):
        """
        Constructs the thesaurus. If we need the ficticious node at the top,
        the only paramenter 'with_root' has to be set to true.
        Left by default said node will not be added
        """
        thesaurus = MeSHThesaurus()
        self.__create_root_nodes(thesaurus)
        #create generic root node to link the different categories
        if with_root:
            self.__add_generic_node(thesaurus)
        #parse the file
        self.__parse_MeSH_file(thesaurus)
        return thesaurus


    #Maps descriptors IDS to tree coordinates.
    def printMappingFile(self, thesaurus,name):
        mapfile = open(name,'w')
        print("Printing mapping file")
        for node_id in thesaurus.get_node_ids():
            node = thesaurus.get_node(node_id)
            #if not node.is_dummy:
            mapfile.write(node_id + "\t" + '\t'.join(node.get_tree_positions()) + "\n")
        print("Done")
        mapfile.close()

    #Maps names to Descriptor IDS
    def printNameIdMapping(self,thesaurus,name):
        mapfile = open(name,'w')
        print("Printing mapping file")
        for node_id in thesaurus.get_node_ids():
            node = thesaurus.get_node(node_id)
            #if not node.is_dummy:
            mapfile.write(node.get_name() + "\t" + '\t'.join(node.get_synonyms()) + '|' + node_id + "\n")
        print("Done")
        mapfile.close()

if __name__ == "__main__":
    #----
    #if len(sys.argv) < 3:
        #print "Usage:"
        #print "python ", sys.argv[0], " descriptors_file annotation_file"
        #sys.exit()
    #----
    print("\t- Loading parser..")
    parser = MeSHParser(sys.argv[1], ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','V','Z'])
    thesaurus = parser.get_thesaurus(True)
    parser.printNameIdMapping(thesaurus, 'full_mapping_synonyms')
    #print "\t- Obtaining annotation"
    #annotation_parser = AnnotationParser(thesaurus, sys.argv[2])
    #annotation = annotation_parser.get_annotations()

