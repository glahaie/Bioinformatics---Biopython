# *- coding:utf-8 -* #

# Pour la question 2: On lit le résultat de l'assemblage par cap3
#, on met les résultats dans un dictionnaire, ensuite on choisit
# 5 contigs au hasard avec un générateur de nombre aléatoire, et
#on blast ces résutats

import random
from Bio.Blast import NCBIWWW

contigs = {}
contig_no = None
contig_seq = ""
contig_size = 0

with open("seq.data.cap.contigs", "r") as f:
    for line in f:
        # on regarde d'abord si c'est un contig ou non
        if line[0] == '>':
            if contig_no == None:
                contig_no = int(line[7:])
            if contig_seq != "":
                contigs.update({contig_no:contig_seq})
                contig_seq = ""
                contig_no = int(line[7:])
                contig_size+= 1
        else :
            contig_seq = contig_seq + line.replace("\n","")
    contigs.update({contig_no:contig_seq})
    contig_size +=1

# Maintenant, on a nos contigs, on en choisit 5 au hasard
random_contig = []

#Je m'assure ici de ne pas avoir de doublon
for i in range(5):
    random_c = random.randint(1, contig_size)
    while random_c in random_contig:
        random_c = random.randint(1,contig_size)
    random_contig.append(random_c)

print random_contig

#On blast maintenant les contigs choisis:
for i in random_contig:
    result_handle = NCBIWWW.qblast("blastn", "nr", contigs[i])

    #on enregistre le résultat
    nom_fichier = "blast_contig_" + str(i) + ".xml"
    save_file = open(nom_fichier, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

print "5 contigs cherchés"



