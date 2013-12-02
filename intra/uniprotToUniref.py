# *- coding:utf-8 *-
#
# uniprotToUniref.py
#
# Par Guillaume Lahaie
# LAHG0407707
#
# Programme qui vérifie si le # obtenu correspond au uniref

from bioservices import UniProt
import sys
import os
from BeautifulSoup import BeautifulSoup

u = UniProt()
with open("uniref_mapping-test.txt", "w") as r:
    with open("uniref_mapping.txt", "r") as f:
        for line in f:
            temp = line.split()
            print("Traitement du contig " + temp[0])
            if temp[1] != "None":
                u.mapping(fr='ACC+ID', to='NF100', query=temp[1])
                res = u.search(temp[1], format='xml', limit=10)
                if res is '':
                    r.write(temp[0] + "\tNone\n")
                else:
                    xml = BeautifulSoup(res)
                    r.write(temp[0] + "\t" + xml.accession.contents[0].strip() + "\n")
            else:
                r.write(temp[0] + "\tNone\n")
