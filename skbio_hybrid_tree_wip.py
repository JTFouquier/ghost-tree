#early draft/pseudocode
#Dictionary improvement over messy code
import os,re,sys,skbio

fin_repset = "rep_set_small.fna"
fin_repset = open(fin_repset,"U")

from skbio.parse.sequences import parse_fasta
repsetdic = {}
for label, seq in parse_fasta(fin_repset,ignore_comment=True):
    repsetdic[label] = seq

fin_taxonomy = "99_otu_taxonomy.txt"
fin_taxonomy = open(fin_taxonomy,"U")

taxdic = {}
taxgendic = {}
for line in fin_taxonomy:
    line = line.split("\t")
    accessionID = line[0]
    taxonomyline = line[1]
    genus = taxonomyline.split(";")
    genus = genus[-2].upper()
    if genus[0:3] == "G__":
        genus = genus[3:]
    taxgendic[accessionID] = genus
    taxdic[accessionID] = taxonomyline

fin_taxonomy.close()
fin_repset.close()

repsetIDlist = []
repsetIDlist = repsetdic.keys()

#print repsetIDlist

#print repsetIDlist[1:30]

repgenlist = []
for i in repsetIDlist:
    genus = taxgendic[i]
    if genus not in repgenlist:
        repgenlist.append(genus)

#print repgenlist


generaSeqIDdic = {}
for m in repgenlist:      #phoma, cladosporium... only from repset
    IDnumlist = []
    generaSeqIDdic[m] = IDnumlist

for key in taxgendic:
    if key in repsetIDlist:
        try:
            g=taxgendic[key]
            generaSeqIDdic[g].append(key)
        except:
            continue


##print len(repsetdic)
##print repsetdic
##print len(generaSeqIDdic)
print generaSeqIDdic

    
#print generaSeqIDdic["USTILAGO"]   #prints accession IDS associated with genus
##print generaSeqIDdic["BIPOLARIS"]

from skbio.sequence import BiologicalSequence


###NEEDS MULTIPLE FILE GENERATION  --> 
fout = open("genus.fasta","w")
for genus in generaSeqIDdic:
    print genus
    seqlist = []
    seqlist = generaSeqIDdic[genus]
    print seqlist
    for i in seqlist:
        seq = repsetdic[i]
        t = BiologicalSequence(seq,id=i)
        line = (t.to_fasta(terminal_character=""))
    fout.write(line)
    fout.write("\n")
    
fout.close()
