# *- coding:utf-8 -*
#
# q5_latex.py
#
# Par Guillaume Lahaie
# LAHG04077707
#
# Dernière modification: 2013-12-19
#
# Assemblage des cartographies des gènes
#

from Bio import SeqIO
from BeautifulSoup import BeautifulSoup
import os

PATH_DIAG = "../diag_genes/"

genes = os.listdir(PATH_DIAG)

nuc = {}
genbank = os.listdir("../genbank_genes/")
for gb in genbank:
    with open("../genbank_genes/"+gb, "r") as f:
        for rec in SeqIO.parse(f, "gb"):
            nuc[gb[:-3]] = len(rec.seq)

with open("q5.tex", "w") as tex:
    for gene in genes:
        tex.write("\\begin{minipage}{\\textwidth}\n")
        tex.write("\\centering{\\bf{\\large{Diagramme du gène " + gene[:-4].replace("_","\\_")+"}}}\n\n")
        if gene[:-4] in nuc:
            tex.write("\\centering{Taille de la séquence: " + str(nuc[gene[:-4]]) +" nucléotides}\n\n")
        tex.write("\\includegraphics[width=\\textwidth]{diag_genes/"+gene+"}\n")
        tex.write("\end{minipage}\n\n")
