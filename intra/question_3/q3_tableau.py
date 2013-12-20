# *- coding:utf-8 -*
#
#q3_tableau.py
#
#Par Guillaume Lahaie
#LAHG04077707
#
#Dernière modification: 2013-12-17
#
#Assemblage des informations pour la création du tableau des contigs pour le numéro 2
#Le formattage en sortie sera celui d'un tableau latex
from Bio import SeqIO
from BeautifulSoup import BeautifulSoup

PATH_GB = "../genbank/"
PATH_BLAST_EBI = "../blastEBI/"
NAME_BLAST_EBI = "blast_result"



#Obtenir la séquence protéique qui a blasté pour le résultat gardé
def obtenirSequenceProt(contig, uniref):
    with open(PATH_BLAST_EBI+NAME_BLAST_EBI+str(contig)+".xml", "r") as xml:
        blast_result = BeautifulSoup(xml)
        temp = blast_result.html.body.ebiapplicationresult.sequencesimilaritysearchresult.hits
        temp2 = temp.findAll('hit')
        for hit in temp2:
            if hit['id'] == uniref:
#On trouve la séquence
                return hit.find('matchseq').contents[0].strip(" \n\t\r")


unirefEBI = {}
with open("../question_2/resultatEBI.txt", "r") as EBI:
    for line in EBI:
        tempEBI = line.split()
        contig_ebi = int(tempEBI[0])
        uniref = tempEBI[1].strip(" \n\t\r")
        seqProt = obtenirSequenceProt(contig_ebi, uniref)
        unirefEBI[contig_ebi] = {"uniref":uniref, "prot":seqProt}

#Au début on traite les numéros d'accession
contigs = {}
with open("../question_2/resultatNBCI.txt", "r") as result:
    for line in result:
        temp = line.split("-|-")
        contig_no = int(temp[0].strip(" \n\t\r"))
        #On va chercher la description en même temps
        gb = open(PATH_GB+temp[2].strip(" \t\n\r") +".gb", 'r')
        for record in SeqIO.parse(gb, "gb"):
            description = record.description
        contigs.update({contig_no:{"genbank":temp[2].strip(" \t\r\n"), "description": description}})

print ("len(contigs) = " + str(len(contigs)))
print ("len(unirefEBI) = " + str(len(unirefEBI)))

with open("contigs_q3.tex", "w") as f:
    f.write("\\footnotesize{\n\\begin{longtable}{|p{1.5cm}|p{2cm}|p{3cm}|p{9cm}|}\n")
    f.write("\\hline\n\\centering{\\bf{Contig}} &")
    f.write("\\centering{\\bf{Accession}} &")
    f.write("\\centering{\\bf{Uniref}} &")
    f.write("\\centering{\\bf{Séquence protéique}}\\\\\n\\endhead\\hline\n")
    keySorted = sorted(unirefEBI.keys())
    for key in keySorted:
        f.write(str(key) + " & ")
        f.write(contigs[key]['genbank'].replace("_", "\\_") + " & ")
        f.write(unirefEBI[key]['uniref'].replace("_", "\\_") + " & ")
        f.write("\\seqsplit{")
        f.write(unirefEBI[key]['prot'])
        f.write("}\\\\\n\\hline\n")
    f.write("\\end{longtable}\n}\n")


