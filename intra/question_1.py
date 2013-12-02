# *- coding:utf-8 -* #

#Question 1: Distribution de taille et de taux de GC
# des contigs. Pour le moment, je calcule les résultats,
#le formattage viendra plus tard.

#J'obtiens les informations directement du fichier zip
#On pourrait aussi modifier pour le faire directement du
#fichier seq.data...

#Pour une représentation graphique, j'associe la longueur du contig à

import zipfile
import sys
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import Gnuplot

contigs = {}
contigs_distribution = {}
contigs_plot = []
histogramme_taux = {}
histogramme_taille = {}
contig_no = None
contig_seq = ""

g = Gnuplot.Gnuplot()

#Dessiner un histogramme
def dessinerHisto (fichierDonnees, xlabel, ylabel, titre, fichierHisto):
    gp = Gnuplot.Gnuplot()
    gp('set ylabel "'+ylabel+'"')
    gp('set xlabel "'+xlabel+'"')
    gp('set xtic rotate by -45 scale 0 font ",8"')

    data = Gnuplot.File(fichierDonnees, using='2:xticlabels(1)', title=titre)
    gp('set style data histograms')
    gp('set style fill solid 1.0 border -1')

    gp.plot(data)
    gp.hardcopy(fichierHisto,terminal = 'postscript', enhanced=1, color=1) #must come after plot() function


#Enregistre le fichier pour les données d'un histogramme
def saveData(data, nom_fichier, step):
    with open(nom_fichier, "w") as histo:
        for k, v in data.iteritems():
            histo.write(str(k*step) + "-"+ str(k*step+step) + "\t"+str(v)+"\n")

#Enregistre le fichier pour les données de min - max - etc
def saveDataMean(data, nom_fichier, cle):
    with open(nom_fichier, "w") as histo:
        for k, v in data.iteritems():
            histo.write(str(k) + "\t"+str(v[cle])+"\n")


#Ajoute les données d'un contig dans tous les endroits nécessaire
def ajouterContig(numero, sequence):
#On garde un tableau de valeurs pour vérification de la distribution
    tailleSeq= len(sequence)
    no_gc = sequence.count("G")
    no_gc += sequence.count("C")

    taux = no_gc / float(tailleSeq) * 100.0
    contigs.update({numero:{"taille":tailleSeq, "tauxGC": taux}})

#On prépare la distribution des tailles et taux GC
    contigs_plot.append([tailleSeq,taux])

#Pour les histogrammes
    insererHisto(taux, 5, histogramme_taux)
    insererHisto(tailleSeq, 10, histogramme_taille)
    return {"taille":tailleSeq, "tauxGC":taux}

#Insérer les valeurs dans les données de l'histogramme
def insererHisto(valeur, step, histogramme):
    cle = int(valeur) / step
    if cle in histogramme:
        histogramme[cle] += 1
    else:
        histogramme[cle] = 1


#Début de 'main'
if len(sys.argv) < 2:
    print "Utilisation: python " + sys.argv[0] + " <archive zip CAP3>"
    exit(1)

totalNucleotide = 0
totalGC = 0
totalContig = 0
with zipfile.ZipFile(sys.argv[1], 'r') as cap3Zip:
    fichier_contigs = cap3Zip.open('seq.data.cap.contigs', 'r')
#Première ligne est toujours un contig
    line = fichier_contigs.readline()
    contig_no = int(line[7:])
    while True:
        line = fichier_contigs.readline()
        if not line or line[0] is '>':
            result = ajouterContig(contig_no, contig_seq)
            totalNucleotide += result["taille"]
            totalGC += result["tauxGC"]
            if not line:
                totalContig = contig_no
                break
            else:
                contig_seq = ""
                contig_no = int(line[7:])
        else:
            contig_seq = contig_seq + line.replace("\n", "")

saveData(histogramme_taux, "histogramme_taux.dat", 5)
saveData(histogramme_taille, "histogramme_taille.dat", 10)

saveDataMean(contigs, "contigs_taille.dat", "taille")
saveDataMean(contigs, "contigs_taux.dat", "tauxGC")

dessinerHisto("histogramme_taux.dat", "Taux de GC (%)", "Nombre de contigs", "Nombre de contigs", "histogramme_taux.eps")
dessinerHisto("histogramme_taille.dat", "Taille du contig", "Nombre de contigs", "Nombre de contigs", "histogramme_taille.eps")

#On appel directement gnuplot - je n'ai pas réussi à obtenir le même
#résultat directement du module python
os.system('gnuplot taux.gp')
os.system('gnuplot taille.gp')

#On calcule la moyenne de taille et de tauxgc
moyenneTaille = float(totalNucleotide) / totalContig
moyenneGC = float(totalGC) / totalContig

#On enregistre un fichier de vérification, qui peut être utilisé directement dans latex
with open("contigs_taille_gc.tex", "w") as f:
    f.write("\\footnotesize{\n\\begin{longtable}{|p{3cm}|p{3cm}|p{3cm}|p{3cm}|p{3cm}|p{3cm}|}\n")
    f.write("\\hline\nContig & Taille & Taux GC & Contig & Taille & Taux GC\\\\\n\hline\n")
    changeLine = False
    for k, v in contigs.iteritems():
        f.write(str(k) + " & " + str(v["taille"]) + "& %.2f" % v["tauxGC"])
        if changeLine:
            f.write("\\\\\n\hline\n")
            changeLine = False
        else:
            f.write(" & ")
            changeLine = True
    if changeLine:
        f.write(" &  & \\\\\n\\hline")
    f.write("\\end{longtable}\n}\n")
    f.write(str(totalContig) + " contigs, taille moyenne: " + str(moyenneTaille) + " Taux GC moyen: " + str(moyenneGC))
print(str(totalContig) + " contigs, taille moyenne: " + str(moyenneTaille) + " Taux GC moyen: " + str(moyenneGC))
