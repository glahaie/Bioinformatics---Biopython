# *- coding: utf-8 * -
#
# q4_genes.py
#
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 19 décembre 2013
#
# Extraction des balises gènes des fichiers
# uniprot


import os
import sys
import re
from BeautifulSoup import BeautifulSoup

PATH_UP = "../uniprot/"

resultat = dict()
intra = os.listdir(PATH_UP)
noGene = open("nogene.txt", "w")
with open("resultatUniProt.txt", "w") as result:
    for xml in intra:
        with open(PATH_UP+xml, "r") as up:
            parsed = BeautifulSoup(up)
            gene = parsed.html.body.find("gene")
            if gene != None:
#On écrit le nom dans le fichier, avec le numéro de Uniprot
                nom = gene.find('name').contents[0].strip(" \t\n\r")
                result.write(xml[:-4] + "\t")
                result.write(nom + "\n")
            else:
                noGene.write(xml[:-4] + "\n")
noGene.close()
#    for line in ncbi:
        #result = line.split("-|-")
        #contig_no = result[0].strip(" \t\n\r")
        #print "traitement de " + contig_no

        #nom_fichier = PATH_EBI + NOM_EBI + contig_no + ".xml"
        #with open(nom_fichier, "r") as result_ebi:
            #xml = BeautifulSoup(result_ebi)
            #temp = xml.html.body.ebiapplicationresult.sequencesimilaritysearchresult.hits
        #if (int(temp['total']) is 0):
            #print("Aucun résultat")
        #else:
            #temp2 = temp.findAll('hit')
            #i = 0
            #for hit in temp2:
                #if TRITICUM in hit['description']:
                    #if MAX_EVALUE > float(hit.find('expectation').contents[0]):
                        #result_file.write(contig_no + "\t" + hit['id'] +  "\n")
                        #break
#result_file.close()
