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
#                i +=1
                #print(str(hit['number']) + " : " + hit["id"])
                #sys.stdout.write("\tscore: " + str(hit.find('bits').contents[0].strip()))
                #sys.stdout.write("\texpectation:" + str(hit.find('expectation').contents[0].strip()))
                #sys.stdout.write("\n\tidentity: " + str(hit.find('identity').contents[0].strip()))
                #sys.stdout.write("\talignment #: " + str(hit.find('alignment')['number']))
                #sys.stdout.write("\tbits: " + str(hit.find('bits').contents[0].strip()) + "\n")
                #if(i > 10):
                    #break
            #while True:
                #choix = raw_input("\n\nChoisir le hit à retenir:  ")
                #match = re.search(r'\d+', choix)
                #if match is None:
                    #print("ERREUR")
                #else:
                    #break
            #choix = int(choix)
            #for hit in temp2:
                #if int(hit['number']) is choix:
                    #print "choix = " + str(choix)
                    #print "uniref = " + hit['id']
                    #resultat.update({contig_no:hit['id']})
                    #break
    #j += 1

##On écrit les résultats dans un fichier
#with open("resultatEBI.txt", 'w') as r:
    #for k, v in resultat.iteritems():
        #r.write(str(k)+"\t"+str(v)+"\n")


