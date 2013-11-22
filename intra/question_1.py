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
contigs_taille = {}
contigs_distribution = {}
contigs_plot_taille = []
contigs_plot_tauxGC = []
contigs_histo = {}
contigs_gc = {}
contig_no = None
contig_seq = ""
contig_size = 0
dataToPlot = []

def insererHisto(taux):
    cle = int(taux) / 5
    if cle in contigs_histo:
        contigs_histo[cle] += 1
    else:
        contigs_histo[cle] = 1

g = Gnuplot.Gnuplot()

if len(sys.argv) < 2:
    print "Utilisation: python " + sys.argv[0] + " <archive zip CAP3>"
    exit(1)

with zipfile.ZipFile(sys.argv[1], 'r') as cap3Zip:
    fichier_contigs = cap3Zip.open('seq.data.cap.contigs', 'r')
    for line in fichier_contigs:
        # on regarde d'abord si c'est un contig ou non
        if line[0] is '>':
            if contig_no == None:
                contig_no = int(line[7:])
            if contig_seq is not "":
                #On a un contig au complet, on peut le traiter
                contigs_taille.update({contig_no:len(contig_seq)})
                contigs_plot_taille.append(len(contig_seq))
                no_gc = contig_seq.count("G")
                no_gc += contig_seq.count("C")
                contigs_gc.update({contig_no:no_gc})
                taux = no_gc / float(len(contig_seq)) * 100.0
                contigs_plot_tauxGC.append(taux)
                contigs_distribution.update({contig_no:taux})
                insererHisto(taux)
#On prépare pour le prochain contig
                contig_seq = ""
                contig_no = int(line[7:])
                contig_size+= 1
        else :
            contig_seq = contig_seq + line.replace("\n","")

# On traite le dernier contig
    contigs_taille.update({contig_no:len(contig_seq)})
    contigs_plot_taille.append(len(contig_seq))
    no_gc = contig_seq.count("G")
    no_gc += contig_seq.count("C")
    contigs_gc.update({contig_no:no_gc})
    taux = no_gc / float(len(contig_seq)) * 100.0
    contigs_plot_tauxGC.append(taux)
    contigs_distribution.update({contig_no:taux})
    insererHisto(taux)

#On dessine un graphique avec gnuplot
d = Gnuplot.Data(contigs_plot_taille, contigs_plot_tauxGC, title='Relation entre la taille des contigs et leur taux de GC')
g('set data style points')
g.xlabel('Taille des contigs')
g.ylabel('Taux de GC')
g.plot(d)
g.hardcopy('points.ps',enhanced=1,color=1)

#On essai un histogramme
gp = Gnuplot.Gnuplot(persist = 1)
gp('set ylabel "Nombre de contigs"')
gp('set xlabel "Pourcentage de GC"')
gp('set xtic rotate by -45 scale 0 font ",8"')
tics = "("
with open("contigs_histogramme.dat", "w") as histo:
    for k, v in contigs_histo.iteritems():
        histo.write(str(k*5) + "-"+ str(k*5+5) + "\t"+str(v)+"\n")

data = Gnuplot.File("contigs_histogramme.dat", using='2:xticlabels(1)', title="Nombre de contig")
#gp(r'set xtics add ' + tics)
gp('set style data histograms')
gp('set style fill solid 1.0 border -1')

gp.plot(data)
gp.hardcopy('histogramme.eps',terminal = 'postscript', enhanced=1, color=1) #must come after plot() function

#Pour le moment, on imprime
with open("contigs_taille_gc.txt", "w") as f:
    f.write("contig\t\ttaille\tTotal GC\tTaux GC\n")
    for k, v in contigs_taille.iteritems():
        f.write("contig " + str(k) + "\t" + str(v) + "\t" + str(contigs_gc[k]) +  "\t" + str(contigs_distribution[k]) + "\n")
