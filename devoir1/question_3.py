# *- coding:utf-8 -* #

# Script pour le numéro 3 du devoir 1: Cette partie ne fait
#qu'envoyer la requête blast au serveur du NCBI, et ensuite
#enregistre le résultat dans un fichier.

from Bio.Blast import NCBIWWW


# Tout d'abord on ouvre le fichier

with open("pANNE.txt", 'r') as f:
    sequence = f.read()

#On enlève les retour de chariot du fichier
sequence.replace("\n","")

#Maintenant, on fait le blast
result_handle = NCBIWWW.qblast("blastn", "nr", sequence)


#on enregistre le résultat
save_file = open("resultat_blast_3.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()
