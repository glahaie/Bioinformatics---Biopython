# *- coding:utf-8 -*
#
# Test pour iterations sur des les fichiers d'un répertoire
import os

intra = os.listdir("intra")
for fichier in intra:
    with open("intra/"+fichier, "r") as lines:
        for line in lines:
            print(line)

#Ça marche pour itérer sur les fichiers, maintenant il faut donc parser le
#xml, et sortir les infos voulues

