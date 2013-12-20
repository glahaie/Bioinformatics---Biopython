# *- coding:utf-8 -*
#
# q2_meanEValue.py
#
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 17 décembre 2013
#
# Analyse des fichiers gb pour obtenir une valeur
# médiane des E-value selon la taille du contig

import sys
import os
import re
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import math


MEAN_VALUE = 109
BLAST_DIR = "../blastNCBI/"

blasts = os.listdir(BLAST_DIR)

#Dictionnaire contenant les résultats
resultLower = {}
resultHigher = {}

def addResult(hsp, dict):
    eValue = int(math.log(hsp.expect, 10))
    if eValue in dict.keys():
        dict[eValue] += 1
    else:
        dict[eValue] = 1


for fichier in blasts:
    with open(BLAST_DIR+fichier, "r") as blast:
        blast_record = NCBIXML.read(blast)
        for alignment in blast_record.alignments:
            if "Triticum" in alignment.hit_def:
#On regarde la longueur du contig
                if blast_record.query_letters < MEAN_VALUE:
#Pour le moment je prends le premier hsp
                    for hsp in alignment.hsps:
                        addResult(hsp, resultLower)
                        break
                else:
                    for hsp in alignment.hsps:
                        addResult(hsp, resultHigher)
                        break
                break
with open("evalue_lower.txt", "w") as f:
    for k,v in resultLower.iteritems():
        f.write(str(k) + " " + str(v) + "\n")
with open("evalue_higher.txt", "w") as f:
    for k,v in resultHigher.iteritems():
        f.write(str(k) + " " + str(v) + "\n")
