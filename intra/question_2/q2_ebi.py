# *- coding:utf-8 *-
#
# q2_ebi.py
#
# Par Guillaume Lahaie
# LAHG0407707
#
# Program qui appel le service rest de ebi pour blaster
# les contigs sur la base de données uniref100

from bioservices import ncbiblast
import zipfile
import sys
import os

path_fichier = "../blastEBI/"
nom_fichier = "blast_result"


if len(sys.argv) < 2:
    print "Utilisation: python " + sys.argv[0] + " <archive zip CAP3>"
    exit(1)

contig_no = 0
contig_seq = ""
contigs = {}

#On prépare les contigs pour le blast
with zipfile.ZipFile(sys.argv[1], 'r') as cap3Zip:
    fichier_contigs = cap3Zip.open('../seq.data.cap.contigs', 'r')
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


s = ncbiblast.NCBIblast()
#On a les données nécessaires pour faire les blast, on le fait
# et on enregistre le résultat du fichier xml
for k, v in contigs.iteritems():
    jobid = s.run(program="blastx", sequence=v,
        stype="dna", database="uniref100", email="glahaie@gmail.com")
    result = s.getResult(jobid, "xml")

    #on enregistre le résultat
    save_file = open(path_fichier+nom_fichier+str(k)+".xml", "w")
    save_file.write(result.prettify())
    save_file.close()
    print("\n\n\nblast du contig " + str(k) + " terminé\n\n\n")

