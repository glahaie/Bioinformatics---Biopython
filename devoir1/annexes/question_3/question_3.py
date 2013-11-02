# *- coding:utf-8 -* #

# Script pour le numéro 3 du devoir 1: Cette partie ne fait
#qu'envoyer la requête blast au serveur du NCBI, et ensuite
#enregistre le résultat dans un fichier.

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

path_fichier = "annexes/question_3/"
nom_resultat = "blast_fichier"
LEN_THRESH = 100
E_THRESH = 1e-50
# Tout d'abord on ouvre le fichier

sequence = ""
with open(path_fichier+"pANNE.txt", 'r') as f:
    for line in f:
        sequence = sequence + line.strip()

#On enlève les retour de chariot du fichier
i = 1
while len(sequence) > LEN_THRESH:

#Maintenant, on fait le blast
    print "i = " + str(i)
    print "on faite un blast sur la séquence de longeur " + str(len(sequence))
    result_handle = NCBIWWW.qblast("blastn", "nr", sequence, megablast=True)


#on enregistre le résultat
    save_file = open(path_fichier+nom_resultat+str(i)+".xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

    list_start = []
    list_end = []
    sequences = []
#Maintenant on enlève de la séquence les zones identifiées
    with open(path_fichier+nom_resultat+str(i)+".xml", "r") as result:
        blast_record = NCBIXML.read(result)
        alignment = blast_record.alignments[0]
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
#On met à jour la séquence pour enlever ce résultat
                if hsp.expect < E_THRESH:
                    list_start.append(hsp.query_start)
                    list_end.append(hsp.query_end)
                
            break
#On a les points à enlever
#sort sur les listes
        list_start.sort()
        list_end.sort()
        start= -1
        end = -1
        for s_start, s_end in zip(list_start, list_end):
            if end < 0:
                sequences.append(sequence[: s_start-1])
                end = s_start-1
            else:
                end = s_start-1
                sequences.append(sequence[start: end])
            start = s_end -1
        sequences.append(sequence[start:])
        sequence = ""
        for fragment in sequences:
            sequence +=fragment
#Pour vérifier les résutats, j'enregistre la nouvelle séquence dans un fichier
        print "On écrit le reste de la séquence avec i = " + str(i)
        with open(path_fichier+"pANNE"+str(i)+".txt", "w") as f:
            f.write(sequence)

    i +=1
