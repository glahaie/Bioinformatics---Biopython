# -* coding:utf-8 *-
#
# q3_hits.py
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 18 décembre 2013
#
# Créé un fichier contenant les régions de hit pour query et subject
# du résultat d'un blast retenu

import sys
import os
import re
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


MAX_EVALUE = 0.01
PATH_BLAST = "../blastNCBI/"
NOM_BLAST = "blast_contig"

resultat = {}

#Ajoute les query_start, query_end, hit_start, hit_end au résultat
def ajouterHSP(contig, hsp, accession):
#On crée des dictionnaires pour hit et query
    hit = {"hit_start": hsp.sbjct_start, "hit_end":hsp.sbjct_end, "qry_start":hsp.query_start, "qry_end":hsp.query_end}

    if not(contig in resultat):
        resultat[contig] = {}
        resultat[contig]["hit"] = []
    resultat[contig]["hit"].append(hit)
    resultat[contig]["frame"] = hsp.frame
    resultat[contig]["accession"] = accession


with open("../question_2/resultatNBCI.txt", "r") as ncbi:
    for line in ncbi:
        temp = line.split("-|-")
        contig_no = int(temp[0].strip(" \n\t\r"))
        accession = temp[2].strip(" \n\t\r")

        with open(PATH_BLAST+NOM_BLAST+str(contig_no)+".xml", "r") as blast:
            blast_record = NCBIXML.read(blast)
            for alignment in blast_record.alignments:
                if alignment.accession == accession:
#On a le bon hit, on regarde les hsps
                    for hsp in alignment.hsps:
                        if hsp.expect < MAX_EVALUE:
                            ajouterHSP(contig_no, hsp, accession)

#On enregistre dans un fichier
with open("hit_locations.txt", "w") as loc:
    for key, value in resultat.iteritems():
        loc.write(str(key) + "\t")
        loc.write(str(value["accession"]) + "\t")
        loc.write(str(value["frame"]) + "\t")
        for hit in value["hit"]:
            loc.write(str(hit["hit_start"]) + "\t")
            loc.write(str(hit["hit_end"]) + "\t")
            loc.write(str(hit["qry_start"]) + "\t")
            loc.write(str(hit["qry_end"]) + "\t")
        loc.write("\n")
