
# -*- coding: UTF-8 -*-
# Atelier 2-2: accesssion à plusieurs aspects, et alignement global (je crois)


#
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline


#Tout d'abord on accède au fichier genbank requis
Entrez.email = "glahaie@gmail.com"
#handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="FR872717")
#seq_record= SeqIO.read(handle, "gb")
#handle.close()


#On écrit dans un fichier
#SeqIO.write(seq_record, "FR872717.gb", "gb")
#On a le bon fichier, maintenant on liste les gènes

#On a le fichier, on va plutôt l'ouvrir

handle = open("FR872717.gb", "rU")
seq_recordIt = SeqIO.parse(handle, "gb")
for seq in seq_recordIt:
    seq_record = seq
handle.close()

#On l'enregistre en fasta maintenant

print "Liste des gènes:"

#Apparamment les features sont dans un dictionnaire où les valeurs sont des listes,
#donc on fait attention, si on veut la chaine, de déréférencer
for feature in seq_record.features:
    if(feature.type == 'CDS'):
        print "Nom du gène: " + feature.qualifiers['gene'][0]
        print "Séquence codante: " + str(feature.location)
        print "Numéro d'accession: " + feature.qualifiers['protein_id'][0]


#On extrait la séquence du gène L1 pour la transformer en proteine
#on a la sous-séquence de l'étape précédente: 
l1 = seq_record.seq[5770:7276]

#maintenant on transcrit
l1_prot = l1.translate(cds=True, table=1)
print l1_prot
l1_prot_seq = SeqRecord(l1_prot)
l1_prot_seq.id = seq_record.id

SeqIO.write(l1_prot_seq, "FR872717.fa", "fasta")


#Maintenant on obtient le fichier correspondant

handle2 = Entrez.efetch(db="protein", rettype="gb", retmode="text", id="CCB84764.1")
prot_record = SeqIO.read(handle2, "gb")
handle2.close()

#on enregistre le fichier
SeqIO.write(prot_record, "CCB84764.fa", "fasta");


#Maintenant on fait l'alignement
needle_cline = NeedleCommandline(asequence="CCB84764.fa", bsequence="FR872717.fa",
        gapopen=5, gapextend=1, outfile="needle.txt")

#Cela produit le fichier en sortie, bizarre
print needle_cline


