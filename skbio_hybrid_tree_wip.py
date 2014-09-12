#early draft/pseudocode
#Dictionary improvement over messy code
import os,re,sys,skbio

"""
variables:
repsetdic = {"JF77364":"ATCGATCGATCG"}
taxdic = {'AY880934': 'k__Fungi;p__Basidiomycota;...;g__Thelephora;s__Thelephora_terrestris\n'}
taxgendic = {'AY880934': 'THELEPHORA'}

repsetIDlist = ['JN032522', 'JF946064', 'JN038494', 'AB044375']
repgenlist = ['UNIDENTIFIED', 'ALTERNARIA', 'FUSARIUM']

generaSeqIDdic ={'FUSARIUM': ['JF773634'], 'ALTERNARIA': ['JN038494', 'JN038495']}

"""

def MakeGeneraFastas(fin_taxonomy,fin_repset):
    fin_repset = open(fin_repset,"U")
    fin_taxonomy = open(fin_taxonomy,"U")
    from skbio.parse.sequences import parse_fasta
    repsetdic = {}
    for label, seq in parse_fasta(fin_repset,ignore_comment=True):
        repsetdic[label] = seq

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

    fin_taxonomy.close()   #taxonomy file from UNITE (this varies!)
    fin_repset.close()

    repsetIDlist = []
    repsetIDlist = repsetdic.keys()

    repgenlist = []
    for i in repsetIDlist:
        genus = taxgendic[i]
        if genus not in repgenlist:
            repgenlist.append(genus)

    generaSeqIDdic = {}
    for m in repgenlist:
        IDnumlist = []
        generaSeqIDdic[m] = IDnumlist

    for key in taxgendic:
        if key in repsetIDlist:
            try:
                g=taxgendic[key]
                generaSeqIDdic[g].append(key)
            except:
                continue

    from skbio.sequence import BiologicalSequence
    for genus in generaSeqIDdic:
        fout = open("g__"+genus+"_seqs.fasta","w")
        seqlist = []
        seqlist = generaSeqIDdic[genus]
        for i in seqlist:
            seq = repsetdic[i]
            t = BiologicalSequence(seq,id=i)
            line = (t.to_fasta(terminal_character=""))
            fout.write(line)
            fout.write("\n")
        fout.close()
    return repsetdic,taxdic,taxgendic,repsetIDlist,repgenlist,generaSeqIDdic

fin_taxonomy = "99_otu_taxonomy.txt"
fin_repset = "rep_set1.fna"

MakeGeneraFastas(fin_taxonomy,fin_repset)

