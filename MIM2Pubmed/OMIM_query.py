
"""
    Simple OMIM 2 PubMed mapper
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



def handleDOM(referenceList,filename):
    mapping = defaultdict(list)
    for reference in referenceList.getElementsByTagName("reference"):
        currentMapping = handleReference(reference)
        if currentMapping:
            mapping[currentMapping[0]].append(currentMapping[1])
    outfile = open(filename,'a')
    for mim,pmed in mapping.items():
        outfile.write(mim+"\t")
        outfile.write('\t'.join(pmed))
        outfile.write("\n")
    outfile.close()

def handleReference(reference):
    pubmed = reference.getElementsByTagName("pubmedID")
    mimNumber = reference.getElementsByTagName("mimNumber")
    #if the pubmed id was not found, we just skip it.
    try:
        return (mimNumber[0].firstChild.data, pubmed[0].firstChild.data)
    except:
        pass


def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return rc

"""
The API will limit how many entries can be retrieved in a single request.
Entries and clinical synopses are limited to 20 per request if any 'includes' are specified, otherwise there is no limit.
Gene map entries are limited to 100 per request.
We add some throttling of our own, just to be nice and make sure we don't get banned.
"""

def fetchData(phenotype_list, outfile, config_file):
    #read the cnf file.
    parser = ConfigParser.ConfigParser()
    parser.read(config_file)
    try:
        api_key = parser.get('APIconfig','key')
        server = parser.get('APIconfig','server')
        time_limit = parser.get('Throttling', 'time')
        req_number = parser.get('Throttling', 'req_number')
    except:
        print('There is a problem with the configuration file. Please refer to the supplementary material for the appropriate format')
        print('Program will now terminate')
        exit()

    #-----
    #just to have a nice progressbar.
    num_lines = sum(1 for line in open(sys.argv[1]))
    bar = progressbar.ProgressBar(maxval = int(num_lines/int(req_number)) + 1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]).start()
    #------------------------
    barCounter = 0

    with open(phenotype_list, 'r') as infile:

        #read the lines.
        lines = infile.readlines()

        current_request_number = False
        current_line = 0

        #check the current second for the throttling
        start_second =  datetime.datetime.now().second

        while (True):
            current_query = 'http://'+server+'/api/entry/referenceList?'
            try:
                selected_lines = lines[current_line : current_line + int(req_number)]
            except:
                selected_lines = lines[current_line:]
            current_line += int(req_number)

            #construct the parameters
            line = ['mimNumber='+ i.strip() + '&' for i in selected_lines]
            #check if not emtpy (for the last line) and then end.
            if not line:
                break;

            current_query = current_query + ''.join(line) + 'apiKey='+api_key
            try:
                response = urllib2.urlopen(current_query)
                dom = xml.dom.minidom.parseString(response.read())
                handleDOM(dom, outfile)
                #just for the bar
                barCounter = barCounter + 1
                bar.update(barCounter)
            except:
                print('The query ' + current_query + ' could not be completed. Full traceback follows')
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)
                exit()

            #throttle. 
            if current_request_number:
                time.sleep(float(req_number))
                start_second = datetime.datetime.now().second
                current_request_number = False
            else:
                current_request_number = True

help_string = """
        Extracts the PubMed identifiers from records in OMIM
        Usage:
        ======
        python OMIM_query.py omim_list_file omim2pubmed_outfile config_file
        \t* omim_list_infile: is a list of the mim numbers to consider.
        \t* omim2pubmed_outfile: output file name
        \t* config_file: is the path for the configuration file with the API details. If left blank ./api_key will be read.
        ---------------------------------------------------------------------------------------------------------------
        """


if __name__ == '__main__':
    
    if len(sys.argv[1:]) < 2:
        print(help_string)
        exit()

    if len(sys.argv) == 2:
        config_file = './api_key'
    else:
        config_file = sys.argv[3]


    fetchData(sys.argv[1], sys.argv[2], config_file)

