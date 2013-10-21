# *- coding:utf-8 -* #


#Parser pour un fichier XML de résultat blast

from Bio.Blast import NCBIXML


nom_fichier = "blast_contig_132.xml"
result_handle = open(nom_fichier)

blast_record = NCBIXML.read(result_handle)

#On choisit une E-VALUE
E_VALUE_THRESH = 0.04

for alignment in blast_record.alignments:
     for hsp in alignment.hsps:
         if hsp.expect < E_VALUE_THRESH:
             print '****Alignment****'
             print 'sequence:', alignment.title
             print 'length:', alignment.length
             print 'e value:', hsp.expect
             print 'total score:',hsp.score
             print hsp.query[0:75] + '...'
             print hsp.match[0:75] + '...'
             print hsp.sbjct[0:75] + '...'

