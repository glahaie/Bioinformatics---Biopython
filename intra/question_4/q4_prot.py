# *- coding:utf-8 *-
#
# q4_prot.py
#
# Par Guillaume Lahaie
# LAHG0407707
#
# Programme qui appel le service rest de uniprot
# pour obtenir des informations sur une séquence représentative

from bioservices import UniProt
import sys
import os

path_fichier = "../uniprot/"

contigs = []

#On obtient l'identifiant uniprot des séquences représentatives
with open("../question_2/resultatEBI.txt", "r") as ebi:
    for line in ebi:
        temp = line.split("\t")
        uniprot = temp[1][10:].strip(" \n\t\r")
        contigs.append(uniprot)

#Maintenant on fait les recherches
u = UniProt()
count = 0
for contig in contigs:
    nom_fichier = path_fichier+contig + ".xml"
    if os.path.isfile(nom_fichier):
        print "Fichier déjà existant: " + contig
    else:
        result = u.searchUniProtId(contig)
        with open(path_fichier+contig+".xml", "w") as uni:
            uni.write(result.prettify())
    count += 1
    print "count = " + str(count) + "contig =  " + contig

print "Nombre de contigs traités: " + str(count)

##Première ligne est toujours un contig
    #line = fichier_contigs.readline()
    #contig_no = int(line[7:])
    #while True:
        #line = fichier_contigs.readline()
        #if not line or line[0] is '>':
            #contigs.update({contig_no:contig_seq})
            #if not line:
                #totalContig = contig_no
                #break
            #else:
                #contig_seq = ""
                #contig_no = int(line[7:])
        #else:
            #contig_seq = contig_seq + line.replace("\n", "")


#s = ncbiblast.NCBIblast()
##On a les données nécessaires pour faire les blast, on le fait
## et on enregistre le résultat du fichier xml
#for k, v in contigs.iteritems():
    #jobid = s.run(program="blastx", sequence=v,
        #stype="dna", database="uniref100", email="glahaie@gmail.com")
    #result = s.getResult(jobid, "xml")

    ##on enregistre le résultat
    #save_file = open(path_fichier+nom_fichier+str(k)+".xml", "w")
    #save_file.write(result.prettify())
    #save_file.close()
    #print("\n\n\nblast du contig " + str(k) + " terminé\n\n\n")

