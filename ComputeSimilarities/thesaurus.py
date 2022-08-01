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
from collections import deque


class ThesaurusNode(object):
    """
    A thesaurus node is a descriptor with the following information:
    a unique `identifier` (alphanumeric), a `name` (alphanumeric), and
    a set of hierarchical relations: a (possibly empty) set of parents
    and a (possibly empty) set of children.
    """
    def __init__(self, name, identifier):
        self.name = name
        self.identifier = identifier
        self.synonyms = list()
        self.parents = set()
        self.children = set()


    def set_synonyms(self, new_synonyms):
        self.synonyms = new_synonyms
    
    def get_synonyms(self):
        return self.synonyms

    def get_name(self):
        return self.name

    def get_identifier(self):
        return self.identifier

    def get_parents(self):
        return self.parents

    def get_children(self):
        return self.children

    def get_ancestors(self):
        queue = deque([self])
        ancestors = set()
        while queue:
            elem = queue.popleft()
            if elem not in ancestors:
                ancestors.add(elem)
                queue.extend(elem.get_parents())
        return ancestors

    def get_descendants(self):
        queue = deque([self])
        descendants = set()
        while queue:
            elem = queue.popleft()
            if elem not in descendants:
                descendants.add(elem)
                queue.extend(elem.get_children())
        return descendants

    def add_parent(self, parent):
        self.parents.add(parent)

    def add_child(self, child):
        self.children.add(child)


class MeSHThesaurusNode(ThesaurusNode):
    """
    A MeSH thesaurus node, containing only information which is interesting
    for us. Apart from the following, the set of trees the node belongs to,
    the positions in the hierarchy (list of strings), and also the set
    of categories (trees and categories are extracted from the list of
    positions)
    """
    def __init__(self, name, identifier, tree_positions):
        super(MeSHThesaurusNode, self).__init__(name, identifier)
        self.tree_positions = tree_positions
        self.__extract_trees_and_categories()
        self.is_dummy = False

    def get_tree_positions(self):
        return self.tree_positions

    """
    In this case, there are two situations that have to be addressed.
    Firstly, the "tree" is the root of the MeSH ontology,i.e. Anatomy, Geographicals,...
    Secondly, the "category", returns the subtree from each tree, i.e. D01,...
    """
    def __extract_trees_and_categories(self):
        self.trees = set()
        self.categories = set()
        for position in self.tree_positions:
            tree = position.split('.')[0]
            self.trees.add(tree)
            category = tree[0]
            self.categories.add(category)

    def set_dummy(self, dummy):
        self.is_dummy = dummy

    def get_dummy(self):
        return self.is_dummy

    def get_trees(self):
        return self.trees

    def get_categories(self):
        return self.categories

    def get_num_categories(self):
        """
        Gives the number of categories of the descriptor. In
        most of the cases it should be equal to one, but it could
        happen that a descriptor belongs to more than one day
        """
        return len(self.categories)


class Thesaurus(object):
    """
    A thesaurus is a set of thesaurus nodes, indexed by their identifier
    """
    def __init__(self):
        self.node_by_id = dict()

    def add_node(self, node):
        self.node_by_id[node.get_identifier()] = node

    def get_node(self, identifier):
        return self.node_by_id[identifier]

    def get_node_ids(self):
        return self.node_by_id.keys()

    def get_nodes(self):
        """
        Gets all the nodes
        """
        return self.node_by_id.values()

    def size(self):
        """
        Returns the number of nodes in the thesaurus
        """
        return len(self.node_by_id)

    def common_ancestors(self, id_1, id_2):
        """
        Retrieves the set of common ancestors given by
        two node ids
        """
        return self.node_by_id[id_1].get_ancestors() & self.node_by_id[id_2].get_ancestors()


class MeSHThesaurus(Thesaurus):
    """
    A MeSH thesaurus is a thesaurus with a set of categories, and
    a set of trees (at least one per category)
    """
    def __init__(self):
        super(MeSHThesaurus, self).__init__()
        self.categories = dict()
        self.trees = dict()
        self.trees_by_category = defaultdict(set)
        self.category_by_tree = dict()
        self.node_by_tree_position = dict()
        self.nodes_by_category = defaultdict(set)

    def add_category(self, category_id, category_name):
        self.categories[category_id] = category_name

    def get_category_ids(self):
        return self.categories.keys()

    def get_nodes_by_category(self, catid):
        return self.nodes_by_category[catid]

    def get_trees_ids(self):
        return self.trees.keys()

    def get_category_name_from_id(self, identifier):
        return self.categories[identifier]

    def get_category_id_from_tree_id(self, identifier):
        return self.category_by_tree[identifier]

    def add_tree(self, tree_id, tree_name):
        self.trees[tree_id] = tree_name
        category_id = (tree_id.split('.')[0])[0]
        self.trees_by_category[category_id].add(tree_id)
        self.category_by_tree[tree_id] = category_id

    def get_node_by_position(self, position):
        return self.node_by_tree_position[position]

    def add_node(self, node):
        super(MeSHThesaurus, self).add_node(node)
        for position in node.get_tree_positions():
            self.node_by_tree_position[position] = node

