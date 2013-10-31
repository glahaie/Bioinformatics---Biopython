#- coding:utf-8 -#

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

gd_diagram = GenomeDiagram.Diagram("Composition du vecteur pANNE.txt")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()
from Bio.SeqFeature import SeqFeature, FeatureLocation

#Essai à la main, à la lecture des fichiers XML

feature = SeqFeature(FeatureLocation(3495, 5341), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)
feature = SeqFeature(FeatureLocation(1, 1058), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)
feature = SeqFeature(FeatureLocation(5798, 6627), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)
feature = SeqFeature(FeatureLocation(5342, 5798), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)
feature = SeqFeature(FeatureLocation(3039, 3494), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)
feature = SeqFeature(FeatureLocation(1742, 2057), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)
feature = SeqFeature(FeatureLocation(2776, 3039), strand = +1)
gd_feature_set.add_feature(feature, name="pHT2", label=True)

gd_diagram.draw(format='linear', pagesize="LETTER", orientation="portrait",fragments=1,
         start=0, end=6227)
gd_diagram.write("GD_labels_default.pdf", "pdf")
