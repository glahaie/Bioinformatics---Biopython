# *-coding:utf-8 *-
# q3_blastx.py
#
#Par Guillaume Lahaie
#LAHG04077707
#
#
#Script qui fait un blastx sur tous les résultats
# de la question 2

import zipfile
import sys
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

path_fichier = "../blastx/"
nom_fichier = "blast_contig"
contigs = {}
contig_seq = ""

if len(sys.argv) < 2:
    print "Utilisation: python " + sys.argv[0] + " <archive zip CAP3>"
    exit(1)

with zipfile.ZipFile(sys.argv[1], 'r') as cap3Zip:
    fichier_contigs = cap3Zip.open('seq.data.cap.contigs', 'r')
#Première ligne est toujours un contig
    line = fichier_contigs.readline()
    contig_no = int(line[7:])
    while True:
        line = fichier_contigs.readline()
        if not line or line[0] is '>':
            contigs.update({contig_no:contig_seq})
            if not line:
                totalContig = contig_no
                break
            else:
                contig_seq = ""
                contig_no = int(line[7:])
        else:
            contig_seq = contig_seq + line.replace("\n", "")

contigs_no = []
#On obtient les numéros de contigs à blaster
with open("hit_locations.txt", "r") as hits:
    for line in hits:
        temp = line.split("\t")
        contigs_no.append(int(temp[0].strip(" \n\t\r")))



#On a les données nécessaires pour faire les blast, on le fait
# et on enregistre le résultat du fichier xml
for k, v in contigs.iteritems():
    if (k in contigs_no) and (k > 71):
        print "traitement du contig " + str(k)
        result_handle = NCBIWWW.qblast("blastx", "nr", v)

        #on enregistre le résultat
        save_file = open(path_fichier+nom_fichier+str(k)+".xml", "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
        print("blast du contig " + str(k) + " terminé")

