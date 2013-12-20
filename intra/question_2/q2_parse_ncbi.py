# *- coding:utf-8 -* #
#
# q2_parse_ncbi.py
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 17 décembre 2013
#
# Parser pour un fichier XML de résultat de blast
# de contigs pour extraire les descriptions et les
# numéros d'accession. Je garde seulement le premier
# correspondant à du blé, ayant une E-value de 0.01 ou moins

import sys
import os
import re
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez


#On choisit une E-VALUE
Entrez.email = "glahaie@gmail.com"
PATH_GB = "../genbank/"
BLAST_DIR = "../blastNCBI/"
E_VALUE_MAX = 0.01

blasts = os.listdir(BLAST_DIR)
resultat = {}

#Ajouter un résultat au fichier de résultat
def ajouterDonnees(alignment, hsp, contig, fichier):
    fichier.write(str(contig) + " -|- " + alignment.title + " -|- " + alignment.accession + "\n")


#Télécharge et enregistre le fichier genbank du numéro d'accession
#Pour le moment ne vérifie pas si le fichier existe déjà
def obtenirFichierGB(accession):
    nom_fichier = PATH_GB+accession+".gb"
    print nom_fichier
    if not os.path.exists(nom_fichier):
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession)
        seq_record = SeqIO.read(handle, "gb")
        handle.close()
        SeqIO.write(seq_record, nom_fichier, "gb")
        print "Fichier enregistré"


with open("resultatNBCI.txt", 'w') as r:
    for fichier in blasts:
        with open(BLAST_DIR+fichier, "r") as f:

            #On trouve le numéro du contig
            tabNo = re.findall(r'\d+', fichier)
            contig_no = int(tabNo[0])

            blast_record = NCBIXML.read(f)
            print("fichier en traitement: " + fichier)
            for alignment in blast_record.alignments:

                if ("Triticum" in alignment.hit_def) and (not alignment.accession[0].isdigit()):
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_MAX:
                            #on garde les données
                            ajouterDonnees(alignment, hsp, contig_no, r)
                            obtenirFichierGB(alignment.accession)
                            break
                    break
