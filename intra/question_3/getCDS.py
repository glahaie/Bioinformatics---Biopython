# *- coding:utf-8 -*
#
# getCDS.py
#
# Par Guillaume Lahaie
# LAHG04077707
#
# Analyse des fichiers gb pour resssortir ceux
# ayant des séquences codantes

from Bio import SeqIO

PATH_GB = "../genbank/"

with open("coding.txt", "w") as result:
    with open("hit_locations.txt", "r") as f:
        for line in f:
            #On va chercher le fichier genbak
            temp = line.split("\t")
            print("traitement du contig " + temp[0])

            handle = open(PATH_GB+temp[1]+".gb", "r")
            for record in SeqIO.parse(handle, "gb"):
                for feature in record.features:
                    if feature.type == "CDS":

#On regarde tous les hits possibles
                        location = feature.location
                        start = location.start
                        end = location.end
                        hit_start = 3
                        hit_end = 4
                        key = "protein_id"
                        while hit_end < len(temp):
                            if (((int(temp[hit_start]) >= start) and (int(temp[hit_start]) <= end)) or
                                ((int(temp[hit_end]) >= start) and (int(temp[hit_end]) <= end))):
#on a un hit
                                result.write(temp[0] + "\t" + temp[1] + "\t")
                                result.write(temp[hit_start] + "\t" + temp[hit_end] + "\t")
                                result.write(temp[hit_start+2] + "\t" + temp[hit_end+2] + "\t")
                                if key in feature.qualifiers:
                                    result.write(feature.qualifiers['protein_id'][0])
                                result.write("\n")

                            hit_start +=4
                            hit_end +=4
