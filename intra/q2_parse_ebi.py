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

resultat = dict()
intra = os.listdir("blastEBI/")
j = 1
for fichier in intra:
    with open("blastEBI/"+fichier, "r") as f:
        print("fichier numéro "+str(j)+": " + fichier)
        tabNo = re.findall(r'\d+', fichier)
        contig_no = int(tabNo[0])
        xml = BeautifulSoup(f)
        temp = xml.html.body.ebiapplicationresult.sequencesimilaritysearchresult.hits
        if (int(temp['total']) is 0):
            print("Aucun résultat")
            resultat.update({contig_no:None})
        else:
            temp2 = temp.findAll('hit')
            i = 0
            for hit in temp2:
                i +=1
                print(str(hit['number']) + " : " + hit["id"])
                sys.stdout.write("\tscore: " + str(hit.find('score').contents[0].strip()))
                sys.stdout.write("\texpectation:" + str(hit.find('expectation').contents[0].strip()))
                sys.stdout.write("\n\tidentity: " + str(hit.find('identity').contents[0].strip()))
                sys.stdout.write("\talignment #: " + str(hit.find('alignment')['number']))
                sys.stdout.write("\tbits: " + str(hit.find('bits').contents[0].strip()) + "\n")
                if(i > 10):
                    break
            while True:
                choix = raw_input("\n\nChoisir le hit à retenir:  ")
                match = re.search(r'\d+', choix)
                if match is None:
                    print("ERREUR")
                else:
                    break
            for hit in temp2:
                if int(hit['number']) is choix:
                    resultat.update({contig_no:hit['id']})
                    break
    j += 1

#On écrit les résultats dans un fichier
with open("resultatEBI.txt", 'w') as r:
    for k, v in resultat.iteritems():
        r.write(str(k)+"\t"+str(v)+"\n")


