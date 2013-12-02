# *- coding:utf-8 *-
#
# gbToUniref.py
#
# Par Guillaume Lahaie
# LAHG0407707
#
# Program qui obtient le uniref correspondant à un no d'accession genbank

from bioservices import UniProt
import sys
import os
from BeautifulSoup import BeautifulSoup

u = UniProt()
with open("uniref_mapping.txt", "w") as r:
    with open("resultatNBCI.txt", "r") as f:
        for line in f:
            temp = line.split()
            print("Traitement du contig " + temp[0])
            u.mapping(fr='EMBL_ID', to='NF100', query=temp[1])
            res = u.search(temp[1], format='xml', limit=10)
            if (int(temp[0]) is 283):
                print res
            if res is '':
                r.write(temp[0] + "\tNone\n")
            else:
                xml = BeautifulSoup(res)
                r.write(temp[0] + "\t" + xml.accession.contents[0].strip() + "\n")

