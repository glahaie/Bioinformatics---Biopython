# *- coding:utf-8 *-

result = []

with open("resultatUniProt.txt", "r") as up:
    for line in up:
        temp = line.split("\t")
        identifiant = temp[0].strip(" \t\n\r")
        if identifiant in result:
            print "Identifiant déjà présent : "+ identifiant
        else:
            result.append(identifiant)
    print result
    print len(result)
