# *- coding:utf-8 -* #
from __future__ import print_function

#Parser pour un fichier XML de résultat blast

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord

#On choisit une E-VALUE
E_VALUE_THRESH = 1e-30
Entrez.email="glahaie@gmail.com"
i = 1
with open("NC000006_blast.xml","r") as result:
    done = []
    blast_record = NCBIXML.read(result)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                if ("Homo sapiens" not in alignment.title) and ("Human" not in alignment.title):
                    if alignment.accession not in done:
                        done.append(alignment.accession)
#On met dans un fichier fasta
                        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=alignment.accession)
                        seq_record = SeqIO.read(handle, "gb")
                        handle.close()
                        SeqIO.write(seq_record, "question4/"+str(i)+"_"+alignment.accession+".gb", "gb")
                        i+=1
                        break
