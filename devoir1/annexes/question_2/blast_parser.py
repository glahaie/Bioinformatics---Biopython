# *- coding:utf-8 -* #


#Parser pour un fichier XML de résultat blast
#Spécifique à la question 2 du devoir. Je sais
#ici que chaque hit a seulement un hsp

import sys
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

#On choisit une E-VALUE
E_VALUE_THRESH = 0.04
Entrez.email = "glahaie@gmail.com"
path = "annexes/question_2/"

path_fichier = path + "blast_contig_"+sys.argv[1] + ".xml"

with open(path_fichier) as fichier:
    blast_record = NCBIXML.read(fichier)
    i = 0
    path_result = path + "contig_"+sys.argv[1]+"/"
    if not os.path.exists(path_result):
        os.makedirs(path_result)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
#On obtient alors le fichier genbank
                handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=alignment.accession, seq_start=hsp.sbjct_start, seq_stop=hsp.sbjct_end)
                seq_record= SeqIO.read(handle, "gb")
                handle.close()
                nom_fichier = path_result + alignment.accession + ".gb"
                SeqIO.write(seq_record, nom_fichier, "gb")

        i += 1
        if i > 10:
            break


