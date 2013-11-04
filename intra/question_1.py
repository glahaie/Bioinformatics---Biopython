# *- coding:utf-8 -* #

#Question 1: Distribution de taille et de taux de GC
# des contigs. Pour le moment, je calcule les résultats,
#le formattage viendra plus tard.

#J'obtiens les informations directement du fichier zip
#On pourrait aussi modifier pour le faire directement du
#fichier seq.data...

import zipfile
import sys
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

contigs = {}
contigs_taille = {}
contigs_gc = {}
contig_no = None
contig_seq = ""
contig_size = 0

if len(sys.argv) < 2:
    print "Utilisation: python " + sys.argv[0] + " <archive zip CAP3>"
    exit(1)

with zipfile.ZipFile(sys.argv[1], 'r') as cap3Zip:
    fichier_contigs = cap3Zip.open('seq.data.cap.contigs', 'r')
    for line in fichier_contigs:
        # on regarde d'abord si c'est un contig ou non
        if line[0] is '>':
            if contig_no == None:
                contig_no = int(line[7:])
            if contig_seq is not "":
                #On a un contig au complet, on peut le traiter
                contigs_taille.update({contig_no:len(contig_seq)})
                no_gc = contig_seq.count("G")
                no_gc += contig_seq.count("C")
                contigs_gc.update({contig_no:no_gc})


#On prépare pour le prochain contig
                contig_seq = ""
                contig_no = int(line[7:])
                contig_size+= 1
        else :
            contig_seq = contig_seq + line.replace("\n","")
# On traite le dernier contig
    contigs_taille.update({contig_no:len(contig_seq)})
    no_gc = contig_seq.count("G")
    no_gc += contig_seq.count("C")
    contigs_gc.update({contig_no:no_gc})

    result_handle = NCBIWWW.qblast("blastn", "nr", contig_seq, megablast=True)
    blast_result = NCBIXML.read(result_handle)
    result_handle.close()
    print blast_result
#On obtient comme preuve de concept le numéro d'accession du premier résultat
    for record in blast_result.alignments:
        for hsp in record.hsps:
            print record.accession
            break
        break

#Pour le moment, on imprime
with open("contigs_taille_gc.txt", "w") as f:
    f.write("contig\ttaille\tTaux GC\nb")
    for k, v in contigs_taille.iteritems():
        f.write("contig " + str(k) + "\t" + str(v) + "\t" + str(contigs_gc[k]) + "\n")
