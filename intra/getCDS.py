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

dir_gb = "genbank/"
with open("resultCDS.txt", "w") as result:
    with open("resultatNBCI.txt", "r") as f:
        for line in f:
            #On va chercher le fichier genbak
            temp = line.split()
            print("traitement du contig " + temp[0])

            handle = open(dir_gb+temp[1]+".gb", "r")
            for record in SeqIO.parse(handle, "gb"):
                for feature in record.features:
                    if feature.type == "CDS":
#On regarde tous les hits possibles
                        location = feature.location
                        start = location.start
                        end = location.end
                        hit_start = 2
                        hit_end = 3
                        key = "protein_id"
                        while hit_end < len(temp):
                            if (((int(temp[hit_start]) >= start) and (int(temp[hit_start]) <= end)) or
                                ((int(temp[hit_end]) >= start) and (int(temp[hit_end]) <= end))):
#on a un hit
                                result.write(temp[0] + "\t" + temp[1] + "\t")
                                result.write(temp[hit_start] + "\t" + temp[hit_end] + "\t")
                                if key in feature.qualifiers:
                                    result.write(feature.qualifiers['protein_id'][0]+"\n")
                                else:
                                    result.write("None\n")

                            hit_start +=4
                            hit_end +=4

