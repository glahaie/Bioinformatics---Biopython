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
from BeautifulSoup import BeautifulSoup

resultat = dict()
intra = os.listdir("blastEBI/")
for fichier in intra:
    with open("blastEBI/"+fichier, "r") as f:
        print("nom de fichier: " + fichier)
        contig_no = int(fichier[12:15].strip('.'))
        print(contig_no)
        xml = BeautifulSoup(f)
        temp = xml.html.body.ebiapplicationresult.sequencesimilaritysearchresult.hits
        temp2 = temp.findAll('hit')
        i = 0
        for hit in temp2:
            i +=1
            print(str(hit['number']) + " : " + hit["id"])
            print("\tscore: " + str(hit.find('score').contents[0].strip()))
            print("\texpectation:" + str(hit.find('expectation').contents[0].strip()))
            print("\tidentity: " + str(hit.find('identity').contents[0].strip()))
            print("\talignment #: " + str(hit.find('alignment')['number']))
            if(i > 10):
                break
    raw_input("\n\nContinuez")



#Ça marche pour itérer sur les fichiers, maintenant il faut donc parser le
#xml, et sortir les infos voulues

