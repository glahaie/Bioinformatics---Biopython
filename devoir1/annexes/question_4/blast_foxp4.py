# *- coding:utf-8 -* #

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

with open("NC_000006_41514164_41570122.gb", "r") as foxp4:
    handle = SeqIO.parse(foxp4, "gb")
    for seq in handle:
#On blast cette séquence
        result_handle= NCBIWWW.qblast('blastn', 'nr', seq.seq, megablast=True)

        save_file=open("NC000006_blast.xml", "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()

print "Résutat du blast sauvegardé dans NC000006_blast.xml"
