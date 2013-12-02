# *- coding:utf-8 -*
#
#q2.tableau.py
#
#Par Guillaume Lahaie
#LAHG04077707
#
#Dernière modification: 2013-12-01
#
#Assemblage des informations pour la création du tableau des contigs pour le numéro 2
#Le formattage en sortie sera celui d'un tableau latex
from Bio import SeqIO


unirefEBI = {}
with open("resultatEBI.txt", "r") as EBI:
    for line in EBI:
        tempEBI = line.split()
        contig_ebi = int(tempEBI[0])
        unirefEBI.update({contig_ebi:tempEBI[1]})

unirefUP = {}
with open("uniref_mapping.txt", "r") as uniref:
    for line in uniref:
        tempUP = line.split()
        contig_up = int(tempUP[0])
        unirefUP.update({contig_up:tempUP[1]})
#Au début on traite les numéros d'accession
contigs = {}
with open("resultatNBCI.txt", "r") as result:
    for line in result:
        temp = line.split()
        contig_no = int(temp[0])
        #On va chercher la description en même temps
        gb = open("genbank/"+temp[1]+".gb", 'r')
        for record in SeqIO.parse(gb, "gb"):
            description = record.description
        contigs.update({contig_no:{"genbank":temp[1], "description": description}})
        if contig_no in unirefEBI:
            contigs[contig_no]["unirefEBI"] = unirefEBI[contig_no]
        if contig_no in unirefUP:
            contigs[contig_no]["unirefUP"] = unirefUP[contig_no]


with open("contigs_q2.tex", "w") as f:
    f.write("\\footnotesize{\n\\begin{longtable}{|p{1.5cm}|p{2cm}|p{6cm}|p{3cm}|p{2cm}|}\n")
    f.write("\\hline\nContig & Accession & Description & Uniref - EBI & Uniref - mapping\\\\\n\hline\n")
    for contig in contigs.iteritems():
        f.write(str(contig[0]) + " & " + contig[1]['genbank'].replace("_", "\\_") + " & " + contig[1]['description'].replace("\\", "\\backslash").replace("_", "\\_") + " & ")
        if "unirefEBI" in contig[1] and contig[1]["unirefEBI"] != "None":
            f.write(contig[1]["unirefEBI"].replace("_", "\\_"))
        f.write(" & ")
        if "unirefUP" in contig[1] and contig[1]["unirefUP"] != "None":
            f.write(contig[1]["unirefUP"].replace("_", "\\_"))
        f.write("\\\\\n\\hline\n")
    f.write("\\end{longtable}\n}\n")


