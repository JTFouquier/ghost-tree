import re
import os
import sys
import math

import skbio

from skbio.parse.sequences import parse_fasta

"""
DIRECTIONS:

This version is no longer depedent on Qiime, but now requires
FastTree (http://www.microbesonline.org/fasttree/)
and MUSCLE (http://www.drive5.com/muscle/)
to be downloaded and installed locally

ALL input files go at the bottom & location of FastTree & Muscle
"""

"""
VARIABLE examples:
silvagenlist = ['UNIDENTIFIED', 'ALTERNARIA', 'FUSARIUM']
repsetdic = {"JF77364":"ATCGATCGATCG"}
taxdic = {'AY880934': 'k__Fungi;p__Basidiomycota;...;g__Thelephora;s__Thelephora_terrestris\n'}
taxgendic = {'AY880934': 'THELEPHORA'}
repsetIDlist = ['JN032522', 'JF946064']
repgenlist = ['UNIDENTIFIED', 'ALTERNARIA', 'FUSARIUM']
generaSeqIDdic ={'FUSARIUM': ['JF773634'], 'ALTERNARIA': ['JN038494', 'JN038495']}
"""

##  Definitions
#########################################

def tic():
    import time
    global start_tic
    start_tic = time.time()
def toc():
    import time
    if 'start_tic' in globals():
        print "Time elapsed:" + str(time.time() - start_tic) + " secs"
    else:
        print "Toc: no tic() found"

def ReduceSilva(silvadb,fixedfastaname):
    global silvadic,silvagenlist
    silvadb = open(silvadb,"U")
    fout = open(fixedfastaname,"w")
    silvadic = {}
    silvagenlist = []
    for label, seq in parse_fasta(silvadb,label_to_name=SilvaLabelChange):                     
        if label == "none":
            continue
        else:
            if label not in silvagenlist:
                silvagenlist.append(label)
                try:
                    silvadic[label] = seq
                    fout.write(">"+label+"\n"+seq+"\n")
                except:
                    continue
    silvadb.close()
    fout.close()
    return silvadic,silvagenlist

def SilvaLabelChange(label):
    try:
        test = re.split(";",label)
    except:
        label = label
        return label
    if re.search("Fungi",label):
        label = re.split(";",label)
        genus = label[-2]
        try:
            genus = re.split(" ",genus)
            genus = genus[0]
        except:
            genus = genus
        if re.search("ceae",genus):
            genus = label[-1]
            try:
                genus = re.split(" ",genus)
                genus = genus[0]
            except:
                genus = genus
        label = genus
        return label
    else:
        label = "none"
        return label
    
def MakeGeneraFastas(fin_taxonomy,fin_repset):
    global repsetdic,taxdic,taxgendic,repsetIDlist,repgenlist,generaSeqIDdic
    fin_repset = open(fin_repset,"U")
    fin_taxonomy = open(fin_taxonomy,"U")
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
        genus = genus[-2]
        if genus[0:3] == "g__":
            genus = genus[3:]
        taxgendic[accessionID] = genus
        taxdic[accessionID] = taxonomyline
    fin_taxonomy.close()
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
                g = taxgendic[key]
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
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if os.path.getsize(file) < 1:
            os.remove(file)
    return repsetdic,taxdic,taxgendic,repsetIDlist,repgenlist,generaSeqIDdic

def AlignToTree():
    #Muscle -> mini ITS alignments
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if file.endswith("_seqs.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+muscledir+" -in "+inputname+" -out " +name+ "_aligned.fasta -quiet -maxiters 2 -diags1")
    #FastTree -> mini ITS trees
    for file in os.listdir(cwd):
        if file.endswith("_aligned.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+ftdir+" -nt -quiet "+inputname+" > " +name+ "_tree.nwk")
 
def InsertITSin18S(backbone,hybridtree):
    #Open the 18S backbone tree file (newick format)
    with open (backbone, "r") as silvafile:
        finaltext = silvafile.read()
    #insert ITS-1 mini trees into 18S backbone
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if file.endswith("_tree.nwk"):
            genusname = str(file)
            genusname = genusname.split("_")
            genusname = genusname[2]
            str(genusname)
            with open (file, "r") as ITS1file:
                ITS1text = ITS1file.read().replace(";","")
                #finaltext = re.sub(r'\bgenusname\b', r'ITS1text', finaltext)
                finaltext = re.sub("\\b"+genusname+"\\b", ITS1text, finaltext)
            #finaltext = finaltext.replace(genusname,ITS1text)
    fout = open(hybridtree,"w")
    finaltext = finaltext.replace(";","")
    finaltext = finaltext.replace("\n","")
    finaltext += ";"
    fout.write(finaltext)
    fout.close()

def CountTotalSeqs():
    cwd = os.getcwd()
    filecount = 0
    global totalcount
    totalcount = 0
    for file in os.listdir(cwd):
        if file.endswith("_seqs.fasta"):
            genusname = str(file)
            genusname = re.split("_",genusname)
            genusname = genusname[2]
            fin = open(file,"U")
            for line in fin:
                if re.search(">",line):
                    filecount +=1
            fin.close()
            totalcount += filecount
            print genusname,"has",filecount,"OTUs"
            filecount = 0
    print "total number of OTUs is ",totalcount
    return totalcount

def CountMissingOTUs():
    global missingOTUs
    filecount = 0
    totalmissing = 0
    missingOTUs = []
    for i in repgenlist:
        if i not in silvagenlist:
            missingOTUs.append(i)
            print i + " from your repset is not found in 18S tree"
    cwd = os.getcwd()
    for i in missingOTUs:
        for file in os.listdir(cwd):
            if file.endswith(i+"_seqs.fasta"):
                genusname = str(file)
                genusname = re.split("_",genusname)
                genusname = genusname[2]
                fin = open(file,"U")
                for line in fin:
                    if re.search(">",line):
                        filecount += 1
                fin.close()
                totalmissing += filecount
                filecount = 0
    print "total number of missing OTUs is ",totalmissing
    accountedfor = float(totalcount-totalmissing)/totalcount*100
    print accountedfor


#########################################
#  Directories & MANDATORY input files
#########################################

#directory of MUSCLE    
muscledir = "/Applications/muscle"
#directory of FASTTREE
ftdir = "/Applications/./FastTree"
#input file
fin_taxonomy = "99_otu_taxonomy.txt"
#input file
fin_repset = "rep_set1.fna"
#input file
silvadb = "SSURef_NR99_115_tax_silva_full_align_trunc.fasta"                          

##  Name of ouput files (optional UNNECESSARY changes)
#########################################

#18S file that has been reduced from large SSU file to include only fungi
fixedfastaname = "SSURef_fixed.fasta"
#Name of 18S backbone phylogenetic tree
backbone = "rep_phylo18Sbackbone.nwk"
#Name of final 18S + ITS hybrid phylogenetic tree
hybridtree = "00_hybridtree.nwk"                             

##  call the functions
#########################################
##tic()
#MakeGeneraFastas(fin_taxonomy,fin_repset)
##print "MakeGeneraFiles: "
##toc()
#ReduceSilva(silvadb,fixedfastaname)
##print "ReduceSilva: "
##toc()
#FixFastaFile(fout,fixedfastaname)
#os.system(""+ftdir+" -nt -quiet "+fixedfastaname+" > "+backbone+"")
#AlignToTree()
#print "AlignToTree: "
#toc()
InsertITSin18S(backbone,hybridtree)
#print "InsertITSin18S: "
#toc()

#These are indeed not even in the SILVA database
##Do not appear to be very abundant fungi (can calculate)
#CountTotalSeqs()
#CountMissingOTUs()

#os.system("filter_tree.py -i 082414_HybridTree.nwk -f rep_set1.fna -o 082414_HybridTree_pruned.nwk") 

######Uncomment these to delete extra files.
#Comment for viewing fasta files, alignment files and tree files.

"""
cwd = os.getcwd()
for file in os.listdir(cwd):
    if file.endswith("_tree.nwk"):
        os.remove(file)
    if file.endswith("_aligned.fasta"):
        os.remove(file)
    if file.endswith("_seqs.fasta"):
        os.remove(file)
        
#os.remove(backbone)
"""
