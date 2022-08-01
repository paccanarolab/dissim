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

#writes the triplet file.
def writeTriplet(file_per_disease,per_disease,sem_sim):
    with open("./localStore/"+ file_per_disease, "w") as f_dis:
        for i in range(per_disease.shape[0]):
            for j in range(i, per_disease.shape[1]):
                if per_disease[i,j] != 0.0:
                    f_dis.write(sem_sim.objects[i] + "\t" + sem_sim.objects[j] + "\t" + str(per_disease[i,j]) + "\n") 

def writeSelectedDescriptor(outfile, values):
    with open("./localStore/"+outfile, "w") as f:
        for key in values:
            for elements in values[key]:
                f.write(str(key) + '\t' +  str(elements[0]) + "\t" +  str(elements[1]) + '\t' + str(elements[2]) + "\n")
