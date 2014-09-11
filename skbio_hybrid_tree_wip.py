# -*- coding: UTF-8 -*-
#early draft/pseudocode
import os,re,sys,skbio
#error no module named skbio.....   
#repset1:

#http://scikit-bio.org/docs/0.2.0/generated/skbio.parse.sequences.parse_fasta.html#skbio.parse.sequences.parse_fasta


from StringIO import StringIO

fin = "rep_set1.fna"
fin = open(fin,"U")

fasta_f = StringIO(fin)

print fasta_f

from skbio.parse.sequences import parse_fasta
for label, seq in parse_fasta(fasta_f,ignore_comment=True):
    print(label,seq)
    break
    

	
	
	





	


