#- coding:utf-8 -#

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
import os
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import zipfile
from Bio.Emboss.Applications import WaterCommandline
from Bio import AlignIO

PATH_GB_GENES = "../genbank_genes/"
PATH_DIAG= "../diag_genes/"
genes = os.listdir(PATH_GB_GENES)
colors = [colors.green,colors.lightgreen, colors.teal, colors.darkgreen,
          colors.seagreen,colors.lawngreen,colors.olivedrab]

#Nouveau diagram, on ajoute la track
def prepareTrack(record, diagram):
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", start=0,
                                                 end=len(record))
    gd_feature_set = gd_track_for_features.new_set()
    return gd_feature_set

def drawExon(feature_set, feature):
    feature = SeqFeature(FeatureLocation(feature.location.start, feature.location.end),
                         strand = +1)
    feature_set.add_feature(feature, name="Exon", label=True, color=colors[0])

def drawContig(contig_no, contig, feature_set):
    feature = SeqFeature(FeatureLocation(contig[0], contig[1]),
                         strand = -1)
    feature_set.add_feature(feature, name="contig "+contig_no, label=True, color="red")
def getStartEnd():
    with open("water.txt", "r") as water:
        firstLine = False
        for line in water:
            if "#" not in line:
                if "seq" in line:
                    if not firstLine:
                        temp = line.split(" ")
                        for word in temp:
                            if word.isdigit():
                                start = int(word)
                                firstLine = True
                                break
                    else:
                        temp = line.split(" ")
                        end = int(temp[-1])
    return (start, end)


def doWater(contig, seq):
    with open("contig.faa", "w") as stuff1:
        stuff1.write(">contig\n")
        stuff1.write(contig)
    with open("seq.faa", "w") as stuff2:
        stuff2.write(">seq\n")
        stuff2.write(str(seq))
    water_cline = WaterCommandline()
    water_cline.asequence="contig.faa"
    water_cline.bsequence="seq.faa"
    water_cline.gapopen=10
    water_cline.gapextend=0.5
    water_cline.outfile="water.txt"
    stdout, stderr = water_cline()
    print(stdout + stderr)
    values = getStartEnd()
    return values

contigs = {}
contig_seq = ""
with zipfile.ZipFile("../Fichiers_CAP3_BLE_PREC.zip", 'r') as cap3Zip:
    fichier_contigs = cap3Zip.open('seq.data.cap.contigs', 'r')
#Première ligne est toujours un contig
    line = fichier_contigs.readline()
    contig_no = int(line[7:])
    while True:
        line = fichier_contigs.readline()
        if not line or line[0] is '>':
            contigs.update({contig_no:contig_seq})
            if not line:
                totalContig = contig_no
                break
            else:
                contig_seq = ""
                contig_no = int(line[7:])
        else:
            contig_seq = contig_seq + line.replace("\n", "")


#Au début, on associe les contigs avec leurs séquences, et aussi les gènes
#avec leurs contigs
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
with open("../question_4/resultatUniProt.txt", "r") as up:
    for line in up:
        temp = line.split("\t")
        uniprot = temp[0].strip(" \t\n\r")
        gene = temp[1].strip(" \t\n\r")
        resultGene[uniprot] = gene


resultGeneContig = {}
for key in resultContigs.keys():
    resultGeneContig[resultGene[key]] = resultContigs[key]

for gene in genes:
    print "traitement de : " + gene
    gd_diagram = GenomeDiagram.Diagram("Gène " + gene[:-3])
    handle = open(PATH_GB_GENES+gene, "r")
    for record in SeqIO.parse(handle, "gb"):
        feature_set = prepareTrack(record, gd_diagram)
        for feature in record.features:
            if feature.type == "exon":
                print "traitement de l'exon :"
                print feature
                print dir(feature)
                print feature.qualifiers["note"]
                drawExon(feature_set, feature)
        #On a dessiné les exons, maintenant on dessine les contigs
        if gene[:-3] in resultGeneContig:
            for contig in resultGeneContig[gene[:-3]]:
                values = doWater(contigs[int(contig)], record.seq)
                drawContig(contig, values, feature_set)
        gd_diagram.draw(format='linear', pagesize=(750, 200), orientation="landscape",fragments=1,
             start=0, end=len(record))
        gd_diagram.write(PATH_DIAG+gene[:-3]+".eps", "eps")
