#- coding:utf-8 -#

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors

gd_diagram = GenomeDiagram.Diagram("Composition du vecteur pANNE.txt")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", start=0, end=6627)
gd_feature_set = gd_track_for_features.new_set()
from Bio.SeqFeature import SeqFeature, FeatureLocation

colors = [colors.green,colors.lightgreen, colors.teal, colors.darkgreen, colors.seagreen,colors.lawngreen,colors.olivedrab]

#Essai à la main, à la lecture des fichiers XML
#blast #1: pHT2
feature = SeqFeature(FeatureLocation(1, 1056), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[0],sigil="ARROW", arrowhead_length=0.5,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(5342, 5798), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[1], sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(3495, 5341), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[2], sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(5798, 6627), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[3], sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(3039, 3494), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[4], sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(1742, 2057), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[5], sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(2776, 3039), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True, color=colors[6], sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)

#blast 2:
feature = SeqFeature(FeatureLocation(1057, 1741), strand = +1)
gd_feature_set.add_feature(feature, name="PGeneClip", label=True, color="blue", sigil="ARROW", arrowhead_length=1)

#blast 3
feature = SeqFeature(FeatureLocation(2058, 2770), strand = +1)
gd_feature_set.add_feature(feature, name="Cloning vector EN.Cherry", label=True, color="red", sigil="ARROW",arrowhead_length=1)


#On essai d'ajouter d'autres track pour représenter la position des blasts
gd_track_for_features = gd_diagram.new_track(1, name="pHT2", start=0, end=4924)
gd_feature_set = gd_track_for_features.new_set()
feature = SeqFeature(FeatureLocation(2246, 4092), strand = None)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[2], sigil="ARROW", arrowhead_length=0.2,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(1, 1058), strand = None)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[0],sigil="ARROW", arrowhead_length=0.2,arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(4093, 4922), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[3], sigil="ARROW", arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(2795, 2339), strand = -1)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[1], sigil="ARROW", arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(2340, 2795), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[4], sigil="ARROW", arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(743, 1058), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[5], sigil="ARROW", arrowshaft_height=0.1)
feature = SeqFeature(FeatureLocation(1984, 2246), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=False, color=colors[6],sigil="ARROW",arrowshaft_height=0.1)

#On essai d'ajouter d'autres track pour représenter la position des blasts
gd_track_for_features = gd_diagram.new_track(1, name="PGeneClip", start=0, end=5267)
gd_feature_set = gd_track_for_features.new_set()
feature = SeqFeature(FeatureLocation(1879, 2563), strand = +1)
gd_feature_set.add_feature(feature, name="PGeneClip", label=False, color="blue", sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)

gd_track_for_features = gd_diagram.new_track(1, name="Cloning vector EN.Cherry", start=0, end=10649)
gd_feature_set = gd_track_for_features.new_set()
feature = SeqFeature(FeatureLocation(7102, 7813), strand = +1)
gd_feature_set.add_feature(feature, name="Cloning Vector EN.Cherry", label=False, color="red", sigil="ARROW", arrowhead_length=1,arrowshaft_height=0.1)

gd_diagram.draw(format='linear', pagesize="LETTER", orientation="portrait",fragments=1,
         start=0, end=10649)
gd_diagram.write("GD_labels_default.eps", "eps")
