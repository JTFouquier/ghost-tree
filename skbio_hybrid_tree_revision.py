import re
import os
import sys
import math

import skbio

from skbio.parse.sequences import parse_fasta

"""
DIRECTIONS:

FastTree (http://www.microbesonline.org/fasttree/)
MUSCLE (http://www.drive5.com/muscle/)


"""

"""
VARIABLE examples:
silvagenlist = ['UNIDENTIFIED', 'ALTERNARIA']
repsetdic = {"JF77364":"ATCGATCGATCG"}
taxdic = {'AY880934': 'k__Fungi;p__Basidiomycota;...;
    g__Thelephora;s__Thelephora_terrestris\n'}
taxgendic = {'AY880934': 'THELEPHORA'}
repsetIDlist = ['JN032522', 'JF946064']
repgenlist = ['UNIDENTIFIED', 'ALTERNARIA', 'FUSARIUM']
generaSeqIDdic ={'FUSARIUM': ['JF773634'], 'ALTERNARIA': ['JN038494', 'JN038495']}
"""

##  Definitions
#########################################

def reduce_silva(silvadb,fixedfastaname):
    """Silva files contain all eukaryotic organisms so this extracts
    only fungal sequences.


    Parameters
    ----------
    silvadic : dict
        A dictionary containing the label (key)
    and sequence (value) from Silva file.

    
    silvagenlist : list
        A list containing all unique fungal genera from Silva file.

    
    

    Returns
    ----------
        silvadic

        silvagenlist

        

    Examples
    ----------





    """

    
    global silvadic,silvagenlist
    silvadb = open(silvadb,"U")
    fout = open(fixedfastaname,"w")
    silvadic = {}
    silvagenlist = []
    for label, seq in parse_fasta(silvadb,label_to_name=silva_label_change):                     
        #change to fungi BOOLEAN
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

def silva_label_change(label):
    """


    Parameters
    ----------


    Returns
    ----------       


    Examples
    ----------



    """

    
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
    """Takes ITS fasta file representative sequences and sorts the
       OTUs/species into their corresponding genus file.


    Parameters
    ----------
    repsetdic : dict
        A dictionary containing the label (key)
        and sequence (value) from ITS representative sequences file.


    repgenlist : list
        A list that contains all unique genera from ITS fasta file.


    taxdic : dict
        A dictionary containing accession ID (key) and the entire
        taxonomy line (value) from the Unite taxonomy file.


    taxgendic: dict
        A dictionary containing accession ID (key) and genus only
        from the Unite taxonomy file. 
    

    repsetIDlist :


    
    


    Returns
    ----------       


    Examples
    ----------



    """
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
        if genus.startswith("g__"):
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
    """
    Parameters
    ----------
    inputname : str
        Entire name of the genus fasta sequence file.
        

    name : str
        Current name only


    
    Returns
    ----------
    files?


    Examples
    ----------
    Many genus fasta files generated, each containing the
    ITS sequence for one or more OTUs/species that fall under
    each genus.  This script turns these fasta files into miniature
    phylogenetic trees that have good alignments within clades.

    Genus fasta file containing two representative seqs:
    

        >AY373854
        AAATTTCCCGGG
        >GQ131879
        AAATTTGGGTTT


    These fasta files will be turned into a multiple sequence
    alignments (MSAs) using the software program
    MUSCLE (http://www.drive5.com/muscle/).
    Multiple sequence files usually contain "-"s, which represent
    a gap in the alignment.

    Genus MSA file:


        >AY373854
        AAATTTCCCGGG---
        >GQ131879
        AAATTT---GGGTTT


    These MSAs will then be turned into phylogenetic trees in
    Newick format using a software program called FastTree.

    Newick files look something like this:


    (EF669954:0.00919,(GU017494:0.03335,HQ702384:0.01090)0.789:0.03099)
    0.999:0.10873)0.950:0.04895)0.808:0.01461)0.983:0.05376)0.820:
    0.00456)0.287:0.00053);

    General overview:

    fasta file -> multiple sequence alignment file -> Newick tree file



    """

    
    #Muscle creates mini ITS alignments
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if file.endswith("_seqs.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+muscledir+" -in "+inputname+" -out " +name+ "_aligned.fasta -quiet -maxiters 2 -diags1")
    #FastTree creates mini ITS trees
    for file in os.listdir(cwd):
        if file.endswith("_aligned.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+ftdir+" -nt -quiet "+inputname+" > " +name+ "_tree.nwk")
 
def InsertITSin18S(backbone,hybridtree):
    """Takes the 18S scaffold tree as a Newick format string and
        replaces the 18S genera locations with the corresponding
        ITS miniature trees, also in Newick format.




    Parameters
    ----------
    finaltext : str
        Starts with 18S Newick data and keeps replacing itself with
        an updated tree with more and more ITS mini tree files.
        Eventually becomes a complete hybrid tree that will be
        written to a file.  


    genusname : str
        Name of the genus from current file in working directory 
        
    
    



    Returns
    ----------       


    Examples
    ----------



    """
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
                ## confirm this is correct?
                finaltext = re.sub("\\b"+genusname+"\\b", ITS1text, finaltext)
    fout = open(hybridtree,"w")
    finaltext = finaltext.replace(";","")
    finaltext = finaltext.replace("\n","")
    finaltext += ";"
    fout.write(finaltext)
    fout.close()

def count_total_seqs():

    """
    Parameters
    ----------
    





    Returns
    ----------       


    Examples
    ----------



    """
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

def count_unplaced_genera():
    """



    Parameters
    ----------
    unplaced_genera : list
        A list that contains the genera (and their OTUs/species) that
        were identified by Unite and which no equivalent genera are
        found in the Silva 18S tree.  Unite has many more fungal
        genera than the Silva 18S database so there
        will be discrepancies.


    filecount : int
        Integer that counts the number of OTUs in one genus file,
        because these OTUs from the genus file
        were unable to be placed into the Silva tree
        due to discrepancies in fungal nomenclature.


    totalmissing : int
        
        

    Returns
    ----------       


    Examples
    ----------



    """
    
    global unplaced_genera
    unplaced_genera = []
    filecount = 0
    totalmissing = 0
    for i in repgenlist:
        if i not in silvagenlist:
            unplaced_genera.append(i)
            print i + " from your repset is not found in 18S tree"
    cwd = os.getcwd()
    for i in unplaced_genera:
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

def remove_accessory_files():
    """ 


    Parameters
    ----------
    cwd : str
        The current working directory
        

    Returns
    ----------       
    

    Examples
    ----------
        Removing or keeping these files depends on whether or not
        the user would like to see each genus' fasta file,
        alignment file, and tree file.  This can be helpful for further
        investigation into specific fungal clades.

    """
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if file.endswith("_tree.nwk"):
            os.remove(file)
        if file.endswith("_aligned.fasta"):
            os.remove(file)
        if file.endswith("_seqs.fasta"):
            os.remove(file)       
    os.remove(backbone)



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
fin_repset = "rep_set1_90percent.fna"
#input file
silvadb = "SSURef_NR99_115_tax_silva_full_align_trunc.fasta"                          

##  Name of ouput files (optional UNNECESSARY changes)
#########################################

#18S file that has been reduced from large SSU file to include only fungi
fixedfastaname = "SSURef_fixed.fasta"
#Name of 18S backbone phylogenetic tree
backbone = "rep_phylo18Sbackbone.nwk"
#Name of final 18S + ITS hybrid phylogenetic tree
hybridtree = "90percentRepSet_hybridtree.nwk"                             

##  call the functions
#########################################
MakeGeneraFastas(fin_taxonomy,fin_repset)
##print "MakeGeneraFiles: "
#reduce_silva(silvadb,fixedfastaname)
os.system(""+ftdir+" -nt -quiet "+fixedfastaname+" > "+backbone+"")
AlignToTree()
InsertITSin18S(backbone,hybridtree)
#print "InsertITSin18S: "


#count_total_seqs()
#count_unplaced_genera()


######Uncomment these to delete extra files.
#Comment for viewing fasta files, alignment files and tree files.



