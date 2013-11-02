# -* coding:utf-8 *-#


#Question 1 - calculer la distribution des nucléotides
# sur la séquence représentant le gène ALS2

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#On lit le fichier

handle = open("NG_008775.gb", "r")
seq_record = SeqIO.parse(handle, 'gb')
for seq in seq_record:
#On parcours la séquence
    dist_a = seq.seq.count("A")
    dist_c = seq.seq.count("C")
    dist_g = seq.seq.count("G")
    dist_t = seq.seq.count("T")
    print "A:  count: " + str(dist_a) + " % = " + \
        str(float(dist_a)/len(seq)*100)
    print "C:  count: " + str(dist_c) + " % = " + \
        str(float(dist_c)/len(seq)*100)
    print "G:  count: " + str(dist_g) + " % = " + \
        str(float(dist_g)/len(seq)*100)
    print "T:  count: " + str(dist_t) + " % = " + \
        str(float(dist_t)/len(seq)*100)
    print "total = " + str(dist_a+dist_c+dist_g+dist_t)
