#!/usr/bin/python
# -*- coding: UTF-8 -*-



#Première partie: on accède au fichier
#Version plus récente, on le fait avec Entrez
#De 5.3.1 - http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc54

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "glahaie@gmail.com"
handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="NM_001123784")
seq_record= SeqIO.read(handle, "gb")
handle.close()

#On écrit dans un fichier
SeqIO.write(seq_record, "NM_001123784.gb", "gb")

#Maintenant on ouvre le fichier et on le parse
for seq_record2 in SeqIO.parse("NM_001123784.gb", "genbank"):
    print "Description: " + seq_record2.description #On trouve le nom de l'espèce ici
    print "Numéro d'accession: " +seq_record2.id #Ou .name
    print "GI: " + seq_record2.name
    print repr(seq_record2.seq)
    print "Longueur de la séquence: %i" % (len(seq_record2)) #Attention, pas comme js!
    print "Mot clés: [%s]" % ', '.join(map(str, seq_record2.annotations['keywords']))
    print "\n\n"
#Pour la taxonomie spécifiquement
    print "taxonomie:"
    print  seq_record2.annotations['taxonomy']

#Version de la séquence:
    print "Version:"
    print seq_record2.annotations['sequence_version']

#Pour le taxon: ne semble pas y avoir de façon facile de le faire, on parse les
#features, et on garde ce qu'on a besoin
    print "taxon:"
    for feature in seq_record2.features:
        if(feature.type == 'source'):
            print feature.qualifiers['db_xref']

#Maintenant, on veut traduire en séquence protéique
#NE PAS OUBLIER: ON DEAL AVEC UNE LISTE, DONC SI ON VEUT LA POSITION
#452 DU TABLEAU, ON ÉCRIT 451
mRNA = seq_record2.seq[451:3496].transcribe()

proteine = mRNA.translate(cds=True)

#Maintenant, en protéine: si on veut seulement une partie de l'ARN, on doit faire une
seq_record3 = SeqRecord(proteine, id='NM_001123784')
seq_record3.description = "Protein"
#Maintenant, on écrit en FASTA
SeqIO.write(seq_record3, "proteine.fasta", "fasta")
