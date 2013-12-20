#- coding:utf-8 -#

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors

gd_diagram = GenomeDiagram.Diagram("Contig")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", start=0, end=1926)
gd_feature_set = gd_track_for_features.new_set()
from Bio.SeqFeature import SeqFeature, FeatureLocation

colors = [colors.green,colors.lightgreen, colors.teal, colors.darkgreen, colors.seagreen,colors.lawngreen,colors.olivedrab]

#Essai à la main, à la lecture des fichiers XML
#blast #1: pHT2
feature = SeqFeature(FeatureLocation(1, 1632), strand = +1)
gd_feature_set.add_feature(feature, name="Exon", label=True, color=colors[0])

#Contig
feature = SeqFeature(FeatureLocation(199, 272), strand = -1)
gd_feature_set.add_feature(feature, name="contig 95", label=True, color="red")

gd_diagram.draw(format='linear', pagesize=(750,200), orientation="landscape",fragments=1,
         start=0, end=1926)
gd_diagram.write("../diag_genes/RFC-1.eps", "eps")
