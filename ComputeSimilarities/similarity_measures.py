from thesaurus import *
from annotation import *
from mesh_parser import MeSHParser
from writeFiles import *
from semsim import *
#--
import numpy as np
import sys

"""
This package contains all semantic similarity measures we have developed.
So far:
    Resnik,
    Lin
    Jiang
    simUI
    simGIC
and Yang's et. al. Random Walk Contribution
"""

class SimUI(SemanticSimilarity):
    """
    Given two diseases i j with extended sets of MeSH terms annotations
    MeSH(i) and MeSH(j) respectively, which include both direct annotations and 
    all their ancestral terms up to the root node, simUI is given by the Jaccard 
    index of the two sets, i.e. the number of terms in the intersection of
    MeSH(i) and MeSH (j) divided by the number of terms in their union
    """
    def __init__(self,thesaurus,annotation):
        self.__thesaurus = thesaurus
        super(SimUI, self).__init__(thesaurus,annotation)

    def semantic_similarity(self, disease1, disease2):
        MeSH1 = set(self.annotation.get_descriptors_per_object(disease1))
        MeSH2 = set(self.annotation.get_descriptors_per_object(disease2))
        value = float(len(MeSH1 & MeSH2)) / float(len(MeSH1 | MeSH2))
        return value

class SimGIC(SemanticSimilarity):
    """
    simGIC is similar to simUI, but it is calculated using a weighted Jaccard index where the
    weight of each element is equal to its information content. simGIC is thus essentially the sum
    of the IC of each term in the intersection of MeSH(i) and MeSH(j), divided by the sum of the IC
    of each term in their unio
    """

    def __init__(self,thesaurus,annotation):
        self.__thesaurus = thesaurus
        super(SimGIC, self).__init__(thesaurus,annotation)

    def semantic_similarity(self, disease1, disease2):
        MeSH1 = set(self.annotation.get_descriptors_per_object(disease1))
        MeSH2 = set(self.annotation.get_descriptors_per_object(disease2))
        intersected = sum([self.annotation.num_annot_per_descriptor(i) for i in (MeSH1 & MeSH2)])
        union = sum([self.annotation.num_annot_per_descriptor(i) for i in (MeSH1 | MeSH2)])
        value = float(intersected) / float(union)
        return value


class Resnik(SemanticSimilarity):
    """
    Resnik (1999) defines semantic similarity between two terms
    c1 and c2 simply as the information content of their most 
    informative common ancestor (i.e. the common ancestor with
    the highest information content)
    """
    def __init__(self,thesaurus,annotation,strategy):
        self.__strategy = strategy
        self.__thesaurus = thesaurus
        super(Resnik, self).__init__(thesaurus,annotation)

    def semantic_similarity(self, id1, id2):

        common_ancestors = self.common_ancestors(id1,id2)
        number_of_annotations  = [self.annotation.num_annot_per_descriptor(i.get_identifier()) for i in common_ancestors]
        min_annot = min(number_of_annotations)
        value = -1.0 * np.log10(float(min_annot) / float(self.num_objects))
        selectedAncestor = common_ancestors[number_of_annotations.index(min_annot)]
        return (selectedAncestor.get_identifier(), value)

    def get_descriptor_indexes(self):
        return super(Resnik, self).get_descriptor_indexes()

    def normalise(self,out):
        return

    def selectionStrategy(self,values):
        if self.__strategy.upper() == "MAX":
            selected_value = np.amax(values)
        if self.__strategy.upper() == "AVG":
            selected_value = np.mean(values)
        if self.__strategy.upper() == "MED":
            selected_value = np.median(values)
        if self.__strategy.upper() == "ALFONSO":
            selected_row = np.amax(values,1)
            selected_col = np.amax(values,0)
            selected_value = np.mean(np.concatenate([selected_row,selected_col]))
        return selected_value

    def get_lowestCommonAncestor(self):
        return super(Resnik,self).get_lowestCommonAncestor()
    
    
class Schlicker(Resnik):

    def __init__(self,thesaurus,annotation,strategy):
        super(Schlicker,self).__init__(thesaurus,annotation,strategy)
        self.__strategy = strategy

    def semantic_similarity(self, id1,id2):
        logp1 = np.log10(float(self.annotation.num_annot_per_descriptor(id1))/float(self.num_objects))
        logp2 = np.log10(float(self.annotation.num_annot_per_descriptor(id2))/float(self.num_objects))
        resnik = super(Schlicker,self).semantic_similarity(id1,id2)
        if logp1 + logp2 == 0:
            lin_similarity = 0
        else:
            lin_similarity = float((-2.0*resnik[1]))/float(logp1+logp2)
        value = lin_similarity * (1 - (np.exp(-resnik[1])))
        return (resnik[0], value)


    #just in case there is some normalisation involved.
    def normalise(self,out):
        return 

    #to define the strategy of selection
    def selectionStrategy(self,values):
        value = super(Schlicker,self).selectionStrategy(values)
        return value


class Lin(Resnik):

    def __init__(self,thesaurus,annotation,strategy):
        super(Lin,self).__init__(thesaurus,annotation,strategy)
        self.__strategy = strategy


    def semantic_similarity(self, id1,id2):
        logp1 = np.log10(float(self.annotation.num_annot_per_descriptor(id1))/float(self.num_objects))
        logp2 = np.log10(float(self.annotation.num_annot_per_descriptor(id2))/float(self.num_objects))
        resnik = super(Lin,self).semantic_similarity(id1,id2)
        if logp1 + logp2 == 0:
            value = 0
        else:
            value = float((-2.0*resnik[1]))/float(logp1+logp2)
        return (resnik[0], value)

    #just in case there is some normalisation involved.
    def normalise(self,out):
        return 

    #to define the strategy of selection
    def selectionStrategy(self,values):
        value = super(Lin,self).selectionStrategy(values)
        return value

class Jiang(Resnik):

    def __init__(self,thesaurus,annotation,strategy):
        super(Jiang,self).__init__(thesaurus,annotation,strategy)
        self.__strategy = strategy

    def semantic_similarity(self,id1,id2):
        logp1 = np.log10(float(self.annotation.num_annot_per_descriptor(id1))/float(self.num_objects))
        logp2 = np.log10(float(self.annotation.num_annot_per_descriptor(id2))/float(self.num_objects))
        resnik = super(Jiang,self).semantic_similarity(id1,id2)
        value = -2.0 * resnik[1] - logp1 - logp2
        return (resnik[0],value)

    def normalise(self,out):
        out = np.subtract(1, np.divide(out,float(np.max(out))))
        return

    def selectionStrategy(self,values):
        value = super(Jiang,self).selectionStrategy(values)
        return value

class ISM(object):

    def __init__(self,thesaurus,annotation,HSM):
        #just stuff we might need
        #ontology and annotations
        self.objects = list(annotation.get_objects())
        self.thesaurus = thesaurus
        self.annotation = annotation
        self.descriptors = list()

        #indices for the objects
        self.__objects_indices = dict(zip(self.objects,range(0,len(self.objects))))
        self.__objects_reverse = dict(zip(range(0,len(self.objects)),self.objects))
        #extract only valid descriptors, i.e. those that have annotations.
        self.descriptors = [i  for i in self.annotation.get_descriptors() if len(self.annotation.get_objects_per_descriptor(i)) > 0]
        #assign a sequential index to all descriptors
        self.__descriptor_indices = dict(zip(self.descriptors,range(0,len(self.descriptors))))
        self.__descriptor_reverse = dict(zip(range(0,len(self.descriptors)),self.descriptors))
        #get those that have zero children, i.e. leaves
        self.__leaves = [self.__descriptor_indices[i] for i in self.descriptors if not self.thesaurus.get_node(i).get_children()]

        #Internals of the ISM algorithm and initialisations
        self.P = np.zeros((len(self.descriptors),len(self.descriptors)))
        self.W = np.matrix(np.identity(len(self.descriptors)))
        self.epsilon = float(0.001)
        self.B = np.zeros((len(self.__leaves), len(self.objects)))
        self.RWC = np.zeros((len(self.objects),len(self.objects)))
        self.HSM = HSM
        self.ISM = np.zeros((len(self.descriptors),len(self.descriptors)))


    def ism(self):
        #initialise matrices
        self.__initialisePmatrix()
        #walk
        self.__walk()
        #initialise B matrix (create A matrix)
        self.__initialiseB()
        ##compute genewise
        self.__genewise()
        ##final combinations
        self.ISM = 0.5 * (self.RWC + self.HSM)
        
    def getISM(self):
        return self.ISM
        

    """Computes the transition matrix. This matrix represents the transition probability from one
    node to another in the DAG."""
    def __initialisePmatrix(self):
        print('Initialising probability matrix..')
        self.P[self.__leaves,self.__leaves] = 1
        #we start in the root and go our way down
        for v in self.descriptors:
            N_c = 0
            N_u = 0
            N_v = len(self.annotation.get_objects_per_descriptor(v))
            N_v_star = self.annotation.get_objects_per_descriptor(v)
            #fetch the children.
            children = [i.get_identifier() for i in self.thesaurus.get_node(v).get_children() if i.get_identifier() in self.descriptors]
            #get the total number of annotations in all the children
            Nc = {}
            for c in children:
                annotations = self.annotation.get_objects_per_descriptor(c)
                Nc[c] = len(annotations)
                N_u = N_u + len(annotations)
                N_v_star = N_v_star - annotations
            #set the individual values
            for c in children:
                N_c = Nc[c]
                self.P[self.__descriptor_indices[c], self.__descriptor_indices[v]] = (1.0-float(len(N_v_star))/float(N_v)) * float(N_c)/float(N_u)
        print('Done!')
    
    def __walk(self):
        print('Walking..')
        self.W_star = self.W
        while True:
            self.W = self.W_star
            self.W_star = self.P * self.W
            matrix_diff = np.linalg.norm(self.W_star - self.W)
            print("\t" + str(matrix_diff) + "/" + str(self.epsilon))
            if matrix_diff <= self.epsilon:
                break
        self.W = self.W_star
        print('Done!')

    def __initialiseB(self):
        print('Computing B..')
        A = self.__computeA()
        #has to do with the indexing of numpy. All columns will be considered. As W is desc x desc, its fine.
        #what does the matrix B represent
        self.B = self.W[self.__leaves] * A
        print('Done!')

    def __computeA(self):
        #transition probability from an annotation to a leaf.
        A = np.zeros((len(self.descriptors), len(self.objects)))
        #for every object.
        for i in self.objects:
            for v in self.annotation.get_descriptors_per_object(i):
                S_i = set([d.get_identifier() for d in self.thesaurus.get_node(v).get_descendants() if len(d.get_children()) == 0])
                A[self.__descriptor_indices[v], self.__objects_indices[i]] = 1.0/float(len(S_i))
        return A

    def __genewise(self):
        #---
        print('Computing RWC..')
        n = len(self.objects)
        sum_col = np.sum(self.B, axis=0).tolist()[0]
        #we should not have zero objects.
        for i in range(n):
            a = np.transpose(self.B[:,i])
            for j in range(i,n):
                b = self.B[:,j]
                combinedSum = a*b
                self.RWC[i][j] = combinedSum / (float(sum_col[i] + sum_col[j]) - combinedSum)
        print("Done!")


if __name__ == "__main__":
    categories = dict([
        ('TWO',['A','C']),
        ('FIVE',['A','C','D','E','G']),
        ('ALL',['A','B','C','D','E','F','G','H','I','J','K','L','M','N','Z'])
        ])
    #----
    if len(sys.argv) < 3:
        print("Usage:")
        print("python ", sys.argv[0], " descriptors_file annotation_file")
        sys.exit()
    #----
    print('Starting utilities..')
    print("\t- Loading parser..")
    parser = MeSHParser(sys.argv[1],categories['ALL'])
    thesaurus = parser.get_thesaurus(False)
    print("\t- Obtaining annotation")
    annotation_parser = AnnotationParser(thesaurus, sys.argv[2])
    #we just use the A category to make it fast
    annotation = annotation_parser.get_annotations('D')

   # for obj in annotation.get_objects():
   #     print obj,
   #     print annotation.get_direct_annotations(obj)
   #     raw_input()
    

    print('Information content per node')
    root = len(list(annotation.get_objects()))
    similarity = Resnik(thesaurus, annotation, 'max')
    with open('./localStore/test_set','w') as f:
        for i in annotation.get_direct_annotations('612460'):
            logp1 = np.log10(float(annotation.num_annot_per_descriptor(i))/float(root))
            value = i + '(' + str(logp1)  + ')\n'
            for j in annotation.get_direct_annotations('610006'):
                common_ancestors = similarity.common_ancestors(i,j)
                number_of_annotations  = [annotation.num_annot_per_descriptor(k.get_identifier()) for k in common_ancestors]
                min_annot = min(number_of_annotations)
                ic = -1.0 * np.log10(float(min_annot) / float(root))
                selectedAncestor = common_ancestors[number_of_annotations.index(min_annot)].get_identifier()
                logp2 = np.log10(float(annotation.num_annot_per_descriptor(j))/float(root))
                value = value + '\t' + j + '('+ str(logp2) +')\t' + selectedAncestor + '(' + str(ic) +')'
                value = value + '\t ' + str(-2.0 * ic/(logp1+logp2)) + '\n'
            f.write(value)
    
