
# *- coding:utf-8 -* #

path_fichier = "annexes/question_3/"
nom_resultat = "blast_fichier"
LEN_THRESH = 100
E_THRESH = 1e-50
# Tout d'abord on ouvre le fichier

sequence = ""
with open(path_fichier+"pANNE.txt", 'r') as f:
    for line in f:
        sequence = sequence + line.strip()


sequence1 = sequence[1056:1742]
sequence2 = sequence[2057:2776]
sequence3 = sequence[3444:3495]

with open(path_fichier+"pANNE1a.txt", "w") as f:
    f.write(sequence1)
with open(path_fichier+"pANNE2a.txt", "w") as f:
    f.write(sequence2)
with open(path_fichier+"pANNE3a.txt", "w") as f:
    f.write(sequence3)

