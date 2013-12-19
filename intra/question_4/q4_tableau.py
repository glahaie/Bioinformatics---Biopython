# *- coding:utf-8 -*
#
# q4_tableau.py
#
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 2013-12-19
#
# Assemblage des informations sur les gènes identifiés à partir des protéines
#

from Bio import SeqIO
from BeautifulSoup import BeautifulSoup

resultContigs = {}
with open("../question_2/resultatEBI.txt", "r") as ebi:
    for line in ebi:
        temp = line.split("\t")
        uniref = temp[1][10:].strip(" \n\t\r")
        if uniref not in resultContigs:
            resultContigs[uniref] = []
        resultContigs[uniref].append(temp[0].strip(" \n\t\r"))

#On prend les données du fichier du nom des gènes
resultGene = {}
with open("resultatUniProt.txt", "r") as up:
    for line in up:
        temp = line.split("\t")
        uniprot = temp[0].strip(" \t\n\r")
        gene = temp[1].strip(" \t\n\r")
        resultGene[uniprot] = gene


with open("contigs_q4.tex", "w") as f:
    f.write("\\footnotesize{\n\\begin{longtable}{|p{5cm}|p{5cm}|p{5cm}|}\n")
    f.write("\\hline\n\\centering{\\bf{Gène}} &")
    f.write("\\centering{\\bf{Référence}} &")
    f.write("\\centering{\\bf{Contigs}} \\\\ \n\\endhead\\hline\n")
    for key in resultGene.keys():
        f.write(str(resultGene[key]).replace("_", "\\_") + " & ")
        if "TRIUR" in resultGene[key]:
            f.write("EnsemblPlants ")
        f.write(" & ")
        for result in resultContigs[key]:
            f.write(result + " ")
        f.write("\\\\\n\\hline\n")
    f.write("\\end{longtable}\n}\n")


