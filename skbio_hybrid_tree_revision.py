import re
import os

import skbio

from skbio.parse.sequences import parse_fasta


cwd = os.getcwd()

"""
DIRECTIONS:

FastTree (http://www.microbesonline.org/fasttree/)
MUSCLE (http://www.drive5.com/muscle/)

"""
def reduce_silva(silvadb):
    """Silva files can contain all eukaryotic organisms so this extracts
       only fungal sequences.  This code will work when someone submits an
       18S database that is not limited to fungi.  It also works for
       18S fasta files that contain only the genus name.

       The 18S input file must be an accurately aligned database such as
       the ones found on Silva.  

    Parameters
    ----------
    fixedfastaname : str
        String that determines that fasta name after only the fungi are
        selected.  
    
    silvadic : dict
        A dictionary containing the label (key) and sequence (value)
        from Silva file.  In try/except block, this tests if label and
        sequence work properly.

    
    silvagenlist : list
        A list containing all unique fungal genera from Silva file.


    label : str
        The name of the DNA sequence corresponding to one DNA sequence.


    seq: str
        The DNA sequence from a fasta file corresponding to a label.

    
    Returns
    ----------
    silvagenlist : list
        A list containing all unique fungal genera from Silva file.
        
    Examples
    ----------
    Below are two smaller examples of fasta sequences from a Silva
    file.  The first one is a fungus and the second one is a plant.
    Only the fungal sequences should be selected (First sequence).


    >AAAA02037088.4053.5842 Eukaryota;Opisthokonta;Nucletmycea;Fungi;
    Dikarya;Ascomycota;Pezizomycotina;Dothideomycetes;
    Dothideomycetidae;Capnodiales;Davidiellaceae;
    Oryza sativa Indica Group
    .......CCU -GGUU----- GA--U-UC-U -G-C-CA--G -UA-G-UCA- ---U--A-U-
    -C-A-AA--- ------G--- AU-U--AA-G --CC-A---- U-G--C---- A-U-G--U-C
    A-UAA--G-C --A------- ---------- ---------- ---------- ----------
    >Eukaryota;Archaeplastida;Chloroplastida;Charophyta;
    Phragmoplastophyta;Streptophyta;Embryophyta;Tracheophyta;
    Spermatophyta;Magnoliophyta;Brassicales;Arabidopsis Arabidopsis	
    .......CCU -GGUU----- GA--U-UC-U -G-C-CA--G -UA-G-UCA- ---U--A-U-
    -C-A-AA--- ------G--- AU-U--AA-G --CC-A---- U-G--C---- A-U-G--U-C
    A-UAA--G-C --A------- ---------- ---------- ---------- ----------


    If the input 18S file has only the genera like below, that is fine
    as well.
    

    >Phoma
    .......CCU -GGUU----- GA--U-UC-U -G-C-CA--G -UA-G-UCA- ---U--A-U-
    -C-A-AA--- ------G--- AU-U--AA-G --CC-A---- U-G--C---- A-U-G--U-C
    A-UAA--G-C --A------- ---------- ---------- ---------- ----------

    
    """    
    global silvagenlist,fixedfastaname
    fixedfastaname = str(silvadb)
    fixedfastaname = re.split("\.",fixedfastaname)
    fixedfastaname = fixedfastaname[0]
    fixedfastaname = fixedfastaname + "_backbone_fungi_for_hybrid_tree.fasta"
    silvadb = open(silvadb,"U")
    fout = open(fixedfastaname,"w")
    silvadic = {}
    silvagenlist = []
    for label,seq in parse_fasta(silvadb,
            label_to_name=silva_label_change):                     
        if label == "Not Fungi":
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
    return silvagenlist

def silva_label_change(label):
    """ Changes the label for the 18S backbone so that it is only the
    genus and nothing more.  


    Parameters
    ----------
    match : str
        String used to identify which input file was submitted (i.e. did
        the input file have genera only or an entire taxonomy line)


    Returns
    ----------
    label : str
        The label or name of each sequence after it is modified.

    Examples
    ----------

    """   
    match = re.split(" ",label)
    try:
        match[1] = ""
    except:
        return label
    if re.search("Fungi",label) and re.search(";",label):
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
        label = "Not Fungi"
        return label
    
def make_genera_fastas(fin_taxonomy,fin_repset):
    """Takes ITS fasta file representative sequences and sorts the
       OTUs/species into their corresponding genus file.  This allows
       OTUs to be compared to other OTUs from the same genus.
       
    Parameters
    ----------
    repsetdic : dict
        A dictionary containing the label (key)
        and sequence (value) from ITS representative sequences file.


    repgenlist : list
        A list that contains all unique genera from ITS fasta file.
        

    taxgendic: dict
        A dictionary containing accession ID (key) and genus only
        from the Unite taxonomy file. ***** not used currently
    

    repsetIDlist : list
        A list that contains all of the IDs from the representative ITS
        sequences.
        

    
    Returns
    ----------       

    Examples
    ----------
    Input is a representative sequence fasta file where each sequence
    corresponds to one representative for all of the OTUs in each
    cluster.  Each sequence has an accession ID that corresponds to
    one sequence in the Unite database.  

    Example of one representative fasta sequence from the input
    fasta file:
    
    
    >>AB015922 Some_comment_ie_sample_location
    CAGAGCCAAGAGATCCGTTGTTGAAAGTTTTTTCAATTCAAGAATAAAACTTAGACTGCAAAG
    ACAACATGAGTTTGGTTTGGGTCTTTGGCGGACACGCTCCAGCCGAAGCCGGTGGGCGGCCGA
    CGCCAGTCCTCACGAACAGCGCCGACGTAGCCCGGCCCGCCAAAGCAACAAGATATAAATCGA
    CACGGGTGGGAGGGTCGACCCAGCACGC


    Example of a taxonomy line:

    AY880934 k__Fungi;p__Basidiomycota;c__Agaricomycetes;
    o__Thelephorales;f__Thelephoraceae;g__Thelephora;
    s__Thelephora_terrestris

    

    This code identifies the genus of all OTUs by looking at the
    accession number from the fasta sequence, then looking at the
    Unite taxonomy file and identifying the genus the sequence
    belongs to. The OTUs then get sorted into genus files that
    have one or more OTUs/species per file.  

    """
    global repgenlist
    fin_repset = open(fin_repset,"U")
    fin_taxonomy = open(fin_taxonomy,"U")
    repsetdic = {}
    for label, seq in parse_fasta(fin_repset,ignore_comment=True):
        repsetdic[label] = seq
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
    for file in os.listdir(cwd):
        if os.path.getsize(file) < 1:
            os.remove(file)
    return repgenlist

def align_to_tree():
    """This function takes each genus file, does a multiple sequence
        alignment with all OTUs from the genus file, and then makes a
        phylogenetic tree in Newick format from the multiple sequence
        alignment.  

    Parameters
    ----------
    inputname : str
        Entire name of the genus fasta sequence file as inputname.
        

    name : str
        Current name used for naming the file

    
    Returns
    ----------
    files?


    Examples
    ----------
    Many genus fasta files are generated, each containing the
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


    for file in os.listdir(cwd):
        if file.endswith("_seqs.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+muscledir+" -in "+inputname+" -out " +name+ "_aligned.fasta -quiet -maxiters 2 -diags1")
    for file in os.listdir(cwd):
        if file.endswith("_aligned.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+ftdir+" -nt -quiet "+inputname+" > " +name+ "_tree.nwk")
 
def insert_ITS_in_18S(backbone,hybridtree):
    """Takes the 18S scaffold tree as a Newick format string and
        replaces the 18S genera locations with the corresponding
        ITS miniature trees, also in Newick format.

    Parameters
    ----------
    finaltext : str
        A string that contains the 18S + ITS hybrid tree information in
        Newick format.


    genusname : str
        Name of the genus from current file in the working directory. 
        

    Returns
    ----------       


    Examples
    ----------
    The fungi from the Silva file are turned into a phylogenetic
    tree in Newick format similar to this with names of all the
    non-redundant genera identified:

    Example of 18S backbone:


    (Phoma:0.00919,(Cladosporium:0.03335,Aspergillus:0.01090)0.789:
    0.03099)0.999:0.10873)0.950:0.04895)0.808:0.01461)0.983:0.05376)
    0.820:0.00456)0.287:0.00053);


    Eexample of ITS genus file with the accession numbers of all
    OTUs/species in the genus file:


    (EF669954:0.00919,(GU017494:0.03335,HQ702384:0.01090)0.789:0.03099)
    0.999:0.10873)0.950:0.04895)0.808:0.01461)0.983:0.05376)0.820:
    0.00456)0.287:0.00053);

    For each genus in the 18S backbone file, the name of each genus
    is replaced with the entire contents of the corresponding genus
    file that contains the ITS miniature phylogenetic tree. So, for
    this example, if there is a genus file for Phoma, Cladosporium or
    Aspergillus, the ITS miniature tree file is inserted in place of
    the genus name in the 18S backbone file.
        

    """
    with open (backbone, "r") as silvafile:
        finaltext = silvafile.read()
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

def count_unplaced_genera():
    """This function counts the OTUs in each genus that were not able to
        be placed onto the 18S scaffold.
    
        Note: Very few OTUs are necessary to get good quality UniFrac 
        results.
        
    Parameters
    ----------
    unplacedgenera: list
        A list of genera within your representative ITS dataset that 
        are not able to be placed into the hybrid tree.
    
    
    totalmissing: int
        Integer that is the total number of OTUs missing from the
        hybrid tree. 
    
    
    fileOTUcount : int
        Integer that is each genus' total number of OTUs
    
    unidentifiedOTUcount : int
        Integer that is the "unidentified genus/group"'s number of OTUs.
        The unidentified group cannot undergo multiple sequence
        alignment, and does not go into the final Silva tree.
        ###### Make sure unidentified does not get into Silva tree
    
    totalcount : int
        Integer that is the total number of OTUs that you have in your
        representative ITS set.
    

    Returns
    ----------       


    Examples
    ----------



    """    
    
    unplacedgenera = []
    fileOTUcount = 0
    totalmissing = 0
    fileOTUcount = 0
    totalcount = 0
    for i in repgenlist:
        if i not in silvagenlist:
            unplacedgenera.append(i)
    for file in os.listdir(cwd):
        if file.endswith("_seqs.fasta"):
            genusname = str(file)
            genusname = re.split("_",genusname)
            genusname = genusname[2]
            fin = open(file,"U")
            for line in fin:
                if re.search(">",line):
                    fileOTUcount +=1
            fin.close()
            totalcount += fileOTUcount
            if genusname == "unidentified":
                unidentifiedOTUcount = fileOTUcount
            if genusname in unplacedgenera:
                totalmissing += fileOTUcount
                print genusname," from repset is not found in 18S tree and has",fileOTUcount,"OTUs"
            fileOTUcount = 0
    print unidentifiedOTUcount," OTUs were unidentified and not placed in tree"
    print "total number of OTUs is ",totalcount
    print "total number of missing OTUs is ",totalmissing
    accountedfor = float(totalcount-totalmissing)/totalcount*100
    print accountedfor

def remove_accessory_files():
    """ Removing or keeping these files depends on whether or not
        the user would like to see each genus' fasta file,
        alignment file, and tree file.  This can be helpful for further
        investigation into specific fungal clades.

    Parameters
    ----------
    cwd : str
        The current working directory
        
    Returns
    ----------       
    
    Examples
    ----------
 
    """
    for file in os.listdir(cwd):
        if file.endswith("_tree.nwk"):
            os.remove(file)
        if file.endswith("_aligned.fasta"):
            os.remove(file)
        if file.endswith("_seqs.fasta"):
            os.remove(file)       
    os.remove(backbone)

#directory of MUSCLE    
muscledir = "/Applications/muscle"
#directory of FASTTREE
ftdir = "/Applications/./FastTree"


fin_taxonomy = "99_otu_taxonomy.txt"
fin_repset = "rep_set1_90percent.fna"
#silvadb = "SSURef_NR99_115_tax_silva_full_align_trunc.fasta"
silvadb = "SSURef_fixed.fasta"

#Name of 18S backbone phylogenetic tree
backbone = "rep_phylo18Sbackbone.nwk"
#Name of final 18S + ITS hybrid phylogenetic tree
hybridtree = "90percentRepSet_hybridtree.nwk"                             


make_genera_fastas(fin_taxonomy,fin_repset)
reduce_silva(silvadb)
os.system(""+ftdir+" -nt -quiet "+fixedfastaname+" > "+backbone+"")
align_to_tree()
insert_ITS_in_18S(backbone,hybridtree)
count_unplaced_genera()
