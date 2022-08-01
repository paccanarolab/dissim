#!/usr/bin/python

"""
    Simple PubMed to MeSH mapper
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

import xml.dom.minidom
import ConfigParser
import urllib2
from itertools import islice
from collections import defaultdict
import datetime,time
import sys,traceback
import progressbar

from Bio import Entrez


class queryPubmed(object):

    def __init__(self, pubmedInputFilename, majorTopic, mapping, pubmedToMeshOutfile, config_file_path):
        self.__pubmedInputFilename = pubmedInputFilename
        self.__majorTopicsOnly = majorTopic
        self.__mapping = mapping
        self.__pubmedToMeshOutfile = pubmedToMeshOutfile
        self.__config_file = config_file_path
        
    def process(self):
        #read the cnf file.
        parser = ConfigParser.ConfigParser()
        parser.read(self.__config_file)
        try:
            config_email = parser.get('EntrezConfig','email')
        except:
            print('There is a problem with the configuration file. Please refer to the supplementary material for the appropriate format')
            print('Program will now terminate')
            exit()
        #-----
        #-------------------------
        #this is just to have a nice progressbar.
        max_simultaneous = 10
        num_lines = sum(1 for line in open(self.__pubmedInputFilename))
        steps = num_lines/max_simultaneous
        if num_lines < max_simultaneous:
            steps = 1
        bar = progressbar.ProgressBar(maxval = steps + 1, widgets=[progressbar.Bar('*', '[', ']'), ' ', progressbar.Percentage()]).start()
        barCounter = 0
        #------------------------
        print("Processing " + str(num_lines) + " records")
        outfile = open(self.__pubmedToMeshOutfile,'a')
        Entrez.email = config_email
        if self.__majorTopicsOnly:
            print('Getting only major topics')

        with open(self.__pubmedInputFilename, 'r') as infile:
            while (True):
                pubmed_ids = [i for i in list(islice(infile, 10))]
                if not pubmed_ids:
                    break;
                if (self.__majorTopicsOnly):
                    #this requires a list
                    allSets = self.getMajorTopics(pubmed_ids)
                else:
                    #this requires a comma separated string.
                    allSets = self.getAllMeSHTerms(','.join(pubmed_ids))

                if len(allSets) > 0:
                    for line in  allSets:
                        if len(line) > 1:
                            outfile.write('\t'.join(line))
                            outfile.write('\n')

                barCounter = barCounter + 1
                bar.update(barCounter)
        bar.finish()
        outfile.close()

    """This function produces lines reading the entire set of pubmed ids"""
    def getAllMeSHTerms(self,pubmed_ids):
        try:
            handle = Entrez.efetch(db="pubmed", id=pubmed_ids, retmode="xml")
        except:
            print('Entrez efetch error. Waiting 5 seconds and retrying...')
            time.sleep(5)
            return self.getAllMeSHTerms(pubmed_ids)
        records = Entrez.read(handle)
        all_mesh_terms = list()
        for record in records:
            try:
                mesh_list = record['MedlineCitation']['MeshHeadingList']
                line = list()
                #add the current pubmed id.
                line.append(str(record['MedlineCitation']['PMID']))
                for mesh_term  in mesh_list:
                    translated = self.translate(str(mesh_term['DescriptorName']))
                    if len(translated) > 0:
                        line.append(translated)
                all_mesh_terms.append(line)
            except:
                pass
            try:
                for rec in records[record]:
                    mesh_list = rec['MedlineCitation']['MeshHeadingList']
                    line = list()
                    #add the current pubmed id.
                    line.append(str(rec['MedlineCitation']['PMID']))
                    for mesh_term  in mesh_list:
                        translated = self.translate(str(mesh_term['DescriptorName']))
                        if len(translated) > 0:
                            line.append(translated)
                    all_mesh_terms.append(line) 
            except:
                continue
                
        return all_mesh_terms

    """
    This function will translate the names provided by the full description of
    the mesh term (i.e. the string) to the unique descriptor. This requires
    an input file with the mapping
    """
    def translate(self,mesh_term):
        translated = str()
        try:
            translated = self.__mapping[mesh_term].strip()
        except KeyError:
            raise
        finally:
            return translated


    """This function produces lines reading only the MajorTopics from pubmed"""
    def getMajorTopics(self,pubmed_ids):
        handle = Entrez.elink(dbfrom="pubmed", db="mesh", id=pubmed_ids)
        records = Entrez.read(handle)
        all_mesh_terms = list()
        for record in records:
            link_sets = record['LinkSetDb']
            line = list()
            line.append(record['IdList'][0])
            for link in link_sets:
                for mesh_headings in link['Link']:
                    line.append(mesh_headings['Id'])
            #write the file.
            all_mesh_terms.append(line)

        return all_mesh_terms

#-----------------------------
def readMappingFile(infile):
    mapping = defaultdict()
    with open(infile,'r') as f:
        for line in f:
            sl = line.strip().split('\t')
            mapping[sl[0]] = sl[1]

    return mapping


help_string = """
        --------------------------------------------------------------------------------------------------------------
        Fetches the MeSH terms from the publications
        Usage:
        ======
        python Pubmed_query.py pubmed_list majorTopicsOnly mappingFile pubmed2mesh_outfile config_file
        \t* pubmed_list: single column file containing the desired PubMed ids.
        \t* majorTopicsOnly: string (Yes/No). Determines whether we get only the major topics (see Supplementary material).
        \t* mappingFile: Double column file, mapping MeSh term names (e.g. Adult) to their unique descriptor identifier (e.g. D000328).
        \t* pubmed2mesh_outfile: writable file where the mappings will be placed.
        \t* config_file: is the path for the configuration file with the API details. If left blank ./entrez_config will be read.
        ---------------------------------------------------------------------------------------------------------------
        """

if __name__ == '__main__':
    
    #just some user output

    if len(sys.argv[1:]) < 3:
        print(help_string)
        exit()

    pubmedInputFile = sys.argv[1]
    pubmed2MeshOutput = sys.argv[4]
    config_file_path = sys.argv[5]

    #read mapping file if necessary
    mapping = defaultdict()
    mappingFile = None
    majorTopicsOnly = True
    if sys.argv[2].upper() == 'NO':
        print('Reading mapping file..')
        mapping = readMappingFile(sys.argv[3])
        majorTopicsOnly = False

    print('Starting to fetch data...')
    #create theobject and process it.
    query = queryPubmed(pubmedInputFile, majorTopicsOnly, mapping, pubmed2MeshOutput, config_file_path)
    query.process()

