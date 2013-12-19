# *- coding:utf-8 -*
#
# q2_parse_ebi.py
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 26 novembre 2013
#
# On itère sur les fichiers de résultat des blasts des contigs
# sur la base de données UNIREF100 pour sortir les hits

import os
import sys
import re
from BeautifulSoup import BeautifulSoup

PATH_EBI = "../blastEBI/"
NOM_EBI = "blast_result"
MAX_EVALUE = 0.01
TRITICUM = "Triticum"

resultat = dict()
intra = os.listdir(PATH_EBI)
result_file = open("resultatEBI.txt", "w")
with open("resultatNBCI.txt", "r") as ncbi:
    for line in ncbi:
        result = line.split("-|-")
        contig_no = result[0].strip(" \t\n\r")
        print "traitement de " + contig_no

        nom_fichier = PATH_EBI + NOM_EBI + contig_no + ".xml"
        with open(nom_fichier, "r") as result_ebi:
            xml = BeautifulSoup(result_ebi)
            temp = xml.html.body.ebiapplicationresult.sequencesimilaritysearchresult.hits
        if (int(temp['total']) is 0):
            print("Aucun résultat")
        else:
            temp2 = temp.findAll('hit')
            i = 0
            for hit in temp2:
                if TRITICUM in hit['description']:
                    if MAX_EVALUE > float(hit.find('expectation').contents[0]):
                        result_file.write(contig_no + "\t" + hit['id'] +  "\n")
                        break
result_file.close()
