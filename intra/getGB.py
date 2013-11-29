# *- coding:utf-8 -*
#
# getGB.py
#
# Par Guillaume Lahaie
# LAHG04077707
#
# Prend les informations du fichier avec le# d'accession, et
# télécharge le fichier genbank


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import os.path

save_path = "genbank/"
Entrez.email = "glahaie@gmail.com"
with open("resultatNBCI.txt", "r") as f:
    for line in f:
        temp = line.split()
        nom_fichier = save_path + temp[1] + ".gb"
        print("Traitement du contig " + str(temp[0]))

        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=temp[1])
        seq_record = SeqIO.read(handle, "gb")
        handle.close()
        SeqIO.write(seq_record, nom_fichier, "gb")

        print("Fichier enregistré")

