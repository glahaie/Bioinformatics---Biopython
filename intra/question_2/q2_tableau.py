# *- coding:utf-8 -*
#
#q2_tableau.py
#
#Par Guillaume Lahaie
#LAHG04077707
#
#Dernière modification: 2013-12-17
#
#Assemblage des informations pour la création du tableau des contigs pour le numéro 2
#Le formattage en sortie sera celui d'un tableau latex
from Bio import SeqIO

PATH_GB = "../genbank/"

unirefEBI = {}
with open("resultatEBI.txt", "r") as EBI:
    for line in EBI:
        tempEBI = line.split()
        contig_ebi = int(tempEBI[0])
        unirefEBI.update({contig_ebi:tempEBI[1]})

#Au début on traite les numéros d'accession
contigs = {}
with open("resultatNBCI.txt", "r") as result:
    for line in result:
        temp = line.split("-|-")
        contig_no = int(temp[0].strip(" \n\t\r"))
        #On va chercher la description en même temps
        gb = open(PATH_GB+temp[2].strip(" \t\n\r") +".gb", 'r')
        for record in SeqIO.parse(gb, "gb"):
            description = record.description
        contigs.update({contig_no:{"genbank":temp[2].strip(" \t\r\n"), "description": description}})
        if contig_no in unirefEBI:
            contigs[contig_no]["unirefEBI"] = unirefEBI[contig_no]

print ("len(contigs) = " + str(len(contigs)))
print ("len(unirefEBI) = " + str(len(unirefEBI)))

with open("contigs_q2.tex", "w") as f:
    f.write("\\footnotesize{\n\\begin{longtable}{|p{1.5cm}|p{2cm}|p{9cm}|p{3cm}|}\n")
    f.write("\\hline\n\\centering{\\bf{Contig}} &")
    f.write("\\centering{\\bf{Accession}} &")
    f.write("\\centering{\\bf{Description}} &")
    f.write("\\centering{\\bf{Uniref - EBI}}\\\\\n\\endhead\\hline\n")
    for contig in contigs.iteritems():
        f.write(str(contig[0]) + " & " + contig[1]['genbank'].replace("_", "\\_") + " & " + contig[1]['description'].replace("\\", "\\backslash").replace("_", "\\_") + " & ")
        if "unirefEBI" in contig[1] and contig[1]["unirefEBI"] != "None":
            f.write(contig[1]["unirefEBI"].replace("_", "\\_"))
        f.write("\\\\\n\\hline\n")
    f.write("\\end{longtable}\n}\n")


