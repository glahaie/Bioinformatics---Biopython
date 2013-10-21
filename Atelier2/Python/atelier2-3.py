# -*- coding: UTF-8 -*-
# Atelier 2-3: Stats sur plusieurs séquences

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline

Entrez.email = "glahaie@gmail.com"
handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="FR872717")
seq_record= SeqIO.read(handle, "gb")
handle.close()

#Donc au début on lit les noms dans le fichier et on les mets dans une liste
liste = open('taxid151340.acc_lst', 'r')

#nouvelle liste
les_ids = []
les_fichiers = []
les_frequences  = {}
for line in liste:
    les_ids.append(line.strip())
liste.close()

#Résultat ok, maintenant on va chercher les infos
Entrez.email = "glahaie@gmail.com"

#ok, on doit le faire un après l'autre
for i in les_ids:
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=i)
    seq_record = SeqIO.read(handle, "gb")
    handle.close()

#On ajoute la clé de ce fichier
    les_frequences[i] = {}
    les_frequences[i]['A'] = seq_record.seq.count('A')
    les_frequences[i]['G'] = seq_record.seq.count('G')
    les_frequences[i]['C'] = seq_record.seq.count('C')
    les_frequences[i]['T'] = seq_record.seq.count('T')
    les_frequences[i]['len'] = len(seq_record)
 
#On calcule les fréquences

for sequence in les_frequences.keys():
    print "fréquences de : " + sequence
    print "A: " + str(float(les_frequences[sequence]['A']) / les_frequences[sequence]['len'])
    print "C: " + str(float(les_frequences[sequence]['C']) / les_frequences[sequence]['len'])
    print "G: " + str(float(les_frequences[sequence]['G']) / les_frequences[sequence]['len'])
    print "T: " + str(float(les_frequences[sequence]['T']) / les_frequences[sequence]['len'])

#Reste les distributions
