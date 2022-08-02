#!/usr/bin/python

import sys
from collections import defaultdict
from itertools import combinations
from rich.progress import track


#MIM:\bGENE_HUMAN\b(UniProt_id),\bGENENAME_HUMAN(UniProtId)...\n
#only three genes in every line, and the line ends with a comma.
def parse_mimtosp(filename, uniprotToGenename, valid_omim):
    mapping = defaultdict(list)
    with open(filename,'r') as f:
        mim = '' 
        for line in f:
            #get rid of the horrible spaces everywere
            line = "".join(line.strip().split())
            pieces = line.strip().split(':')
            if len(pieces) ==  2:
                mim = line.split(':')[0].strip()
                genes = line.split(':')[1].strip()
            else:
                genes = line

            for gene in genes.split(','):
                if not gene:
                    break

                start = str(gene).index('(') + len('(')
                end = str(gene).index(')', start)

                if mim in valid_omim:
                    mapping[mim].append(gene[start:end])
                    mapping[mim].extend(uniprotToGenename[gene[start:end]])
    return mapping


def parse_validOmim(filename):
    valid_omim = dict()
    with open(filename,'r') as f:
        for line in f:
            valid_omim[line.strip()] = 1
    return valid_omim


def parse_uniprotToGenename(filename):
    mapping = defaultdict(list)
    with open(filename,'r') as f:
        for line in f:
            pieces = line.split('\t')
            symbol = pieces[0].strip()
            synonyms = pieces[1]
            uniprot = pieces[2].strip()
            if uniprot == '':
                continue
            mapping[uniprot].append(symbol)
            if len(synonyms) > 0:
               mapping[uniprot].extend([i.strip() for i in synonyms.split(',')])
    return mapping

#=========================================================================

class Interaction(object):
    __columnsHPRD__ = [ "Interactor_1_Gene_symbol", "Interactor_1_HPRD_id", "Interactor_1_RefSeq_id", "Interactor_2_Gene_symbol", "Interactor_2_HPRD_id",
                        "Interactor_2_RefSeq_id", "Experiment_type", "Pubmed_ReferenceList"]
    __columnsHINT__ = ["Interactor_1_HINT_id", "Interactor_2_HINT_id", "Interactor_1_Gene_symbol", "Interactor_2_Gene_symbol", "PubMed_ReferenceList"]

    def __init__(self, source,  *args):
        self.__mapping = defaultdict()
        if len(args) == 1:
            args = args[0].strip().split("\t")
        for (name, value) in zip(getattr(self,"__"+source+"__"),args):
            setattr(self, name, value)

class parser(object):
    def __init__(self, filename, source):
        self.__filename = filename
        self.__source = source

    def mapping(self):
        with open(self.__filename,'r') as mapFile:
            for line in mapFile:
                try:
                    yield Interaction(self.__source, line)
                except TypeError as ex:
                    print("Parse Error. Please check the file")
                    exit()
    
    def __iter__(self):
        return self.mapping()


#=========================================================================

class Interactions(object):
    """columnsHPRD"""
    """columnsHINT"""
    def __init__(self,filename,fileDescription):
        self.__parser = parser(filename, fileDescription)
        self.__ppi = defaultdict(list)
        self.__loadInteractions()

    def __loadInteractions(self):
        for interaction in self.__parser:
            self.__ppi[interaction.Interactor_1_Gene_symbol].append(interaction.Interactor_2_Gene_symbol)
            if interaction.Interactor_1_Gene_symbol != interaction.Interactor_2_Gene_symbol:
                self.__ppi[interaction.Interactor_2_Gene_symbol].append(interaction.Interactor_1_Gene_symbol)

    def getInteractions(self):
        return self.__ppi

#=========================================================================


def produce_benchmark(mimtosp, ppiNetwork, outfilename, useSharedProteins=True):
    allPairs = combinations(mimtosp.keys(),2)
    allPairslength = len(mimtosp.keys())
    val = ((allPairslength)*(allPairslength-1))//2
    with open(outfilename, 'w') as f:
        leftValue = ''
        rightValue = ''
        for diseasePair in track(allPairs, total=val, description="Producing benchmark..."):
            #just to make sure we allways write the smaller value to the left.
            leftValue = min(diseasePair)
            rightValue = max(diseasePair)
            #the diseases share a protein. They need only to share a single protein
            #for the pair to be excluded.
            if len(set(mimtosp[leftValue]) & set(mimtosp[rightValue])):
                #to avoid the cases where proteins are shared.
                if  useSharedProteins:
                    f.write(leftValue +"\t" + rightValue + "\t" + str(1) + "\n")
                else:
                    continue
            #check the interactions
            # we are looking for diseases whose proteins interact with one another.
            else:
                #fetch the interactors of the first disease
                interactors = list()
                for gene in mimtosp[diseasePair[1]]:
                    interactors.extend(ppiNetwork[gene])
                #check if they interact with one of the other proteins.
                if set(interactors) & set(mimtosp[diseasePair[0]]):
                    f.write(leftValue +"\t" + rightValue + "\t" + str(1) + "\n")
                else: 
                    #the reverse case.
                    interactors = list()
                    for gene in mimtosp[diseasePair[0]]:
                        interactors.extend(ppiNetwork[gene])
                    if set(interactors) & set(mimtosp[diseasePair[1]]):
                        f.write(leftValue +"\t" + rightValue + "\t" + str(1) + "\n")

help_string = """
--------------------------------------------------------------------------------------------------------------
This script will produce a matrix (in triplet format) of
MIM numbers that share interactors in the ppi network.
Since there are very many data sets, we provided a mechanism
to generate the data from these various datasets.
Input:
        Produces PPIN
        Usage:
        ======
        python MIM2gene.py uniprot_to_genename mim2sp ppi_file phenotypes outfile
        Where:
        \t*HGNC name to UniProt mapping file.
        \t*mimtosp file (http://www.uniprot.org/docs/mimtosp.txt)
        \t*protein-protein interaction file
        \t*valid_omim file of accepted omim numbers, e.g. phenotype list file.
        \t*output_file_name: A path for the output file.
--------------------------------------------------------------------------------------------------------------
"""

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print(help_string)
        exit()
    uniprotToGenename = parse_uniprotToGenename(sys.argv[1])
    valid_omim = parse_validOmim(sys.argv[4])
    mimtosp = parse_mimtosp(sys.argv[2],uniprotToGenename,valid_omim)
    hprd = Interactions(sys.argv[3], 'columnsHPRD')
    #produce benchmark
    produce_benchmark(mimtosp, hprd.getInteractions(), sys.argv[5], False)
