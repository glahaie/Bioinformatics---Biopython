# *- coding:utf-8 *-
#
# gbToUniref.py
#
# Par Guillaume Lahaie
# LAHG0407707
#
# Dernière modification: 17 décembre 2013
#
# Program qui obtient le uniref correspondant à un no d'accession genbank

from bioservices import UniProt
import sys
import os
from BeautifulSoup import BeautifulSoup

UNIREF_PATH = "../uniref/"

u = UniProt()
with open("uniref_mapping.txt", "w") as r:
    with open("resultatNBCI.txt", "r") as f:
        for line in f:
            temp = line.split("-|-")
            print("Traitement du contig " + temp[0])
            accession = temp[2].strip(" \t\n\r")
            u.mapping(fr='EMBL_ID', to='NF100', query=accession)
            res = u.search(accession, format='xml', limit=10)
            if res is '':
                r.write(temp[0] + "\tNone\n")
                print "aucun résultat pour ce contig"
            else:
                contig = temp[0].strip(" \t\n\r")
                with open(UNIREF_PATH+"result"+contig+".xml", "w") as xml:
                    xml.write(res)
                #xml = BeautifulSoup(res)
                r.write(contig + "\t Result\n")

