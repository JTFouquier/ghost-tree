import re,os,sys,skbio

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

def FixFastaFile(fastafile,fixedfastaname):
    fastafile = open(fastafile,"U")
    fixedfastaname = open(fixedfastaname,"w")
    for line in fastafile:
        if re.match(r'^\s*$', line):
            continue
        if re.search(">",line):
            line = line.strip()
            fixedfastaname.write(line)
            fixedfastaname.write("\n")
        else:
            line = re.sub(" ","",line)
            line = re.sub("\.","-",line)
            line.strip()
            fixedfastaname.write(line)
    fastafile.close()
    fixedfastaname.close()


"""
##definition in SKBIO causes encode/decode error...??
def ReduceSilva(silvadb):  #non redundant, fungi only in my list? faster
    count = 0
    from skbio.parse.sequences import parse_fasta
    FixFastaFile(silvadb,"silva_med_fixed.txt")
    silvadb2 = open("silva_med_fixed.txt","U")
    silvadic = {}
    for label, seq in parse_fasta(silvadb2,label_to_name=SilvaLabelChange):                     
        if label == "none":
            print "NO fungi"
        else:
            try:
                silvadic[label] = seq
                print "YES fungi"
            except:
                print label
        count+=1
        if count == 500:
            break
    silvadb2.close()
    return silvadic

def SilvaLabelChange(label):
    if re.search("Fungi",label):
        label = re.split(";",label)
        genus = label[-2]
        if re.search("ceae",genus):
            genus = label[-1]
        label = genus
    else:
        label = "none"
    return label
"""
                    
def ReduceSilva(silvadb,fout):
    silvadb = open(silvadb,"U")
    fout1 = "intermediate1.fasta"
    fout1 = open(fout1,"w")
    fungi = 0
    for line in silvadb:
        if re.search(">",line):
            line = re.sub(" ","_",line)
        line = re.sub(" ","",line)
        line = line.strip()
        if re.search(">",line) and re.search("Fungi",line):
            fout1.write("\n")
            fout1.write(line)
            fout1.write("\n")
            fungi = 1
        if fungi == 1 and not re.search("Fungi",line) and not re.search(">",line):
            fout1.write(line)
        if re.search(">",line) and not re.search("Fungi",line):
            fungi = 0
    silvadb.close()
    fout1.close()
    fin2 = open("intermediate1.fasta","U") #intermediate1 becomes input
    fout2 = "intermediate2.fasta"
    fout2 = open(fout2,"w")
    global silvagenlist
    silvagenlist = []
    alignmentseq = ""
    genusname = ""
    for line in fin2:
        if re.search("-------",line):
            alignmentseq = line
        if re.search(">",line):
            #line = re.sub("_Eukaryota","",line)
            line = line.split(";")
            genusorfamily = line[-2]
            genusorfamily = genusorfamily[-4:]
            if genusorfamily == "ceae":
                genusname = line[-1]
                genusname = re.split("_",genusname)
                genusname = genusname[0]
            if genusorfamily != "ceae":
                genusname = line[-2]
            if genusname not in silvagenlist:
                silvagenlist.append(genusname)
                fout2.write(alignmentseq)
                fout2.write("\n")
                genusname = ">" + genusname
                fout2.write(genusname)
                fout2.write("\n")
    fout2.write(alignmentseq)               
    fin2.close()
    fout2name = fout2.name
    fout2.close()
    os.remove("intermediate1.fasta")
    FixFastaFile(fout2name,fout)
    os.remove("intermediate2.fasta")
    return silvagenlist
    
def Warnings():
    for i in repgenlist:
        if i not in silvagenlist:
            print i + " from your repset is not found in 18S tree"

def MakeGeneraFastas(fin_taxonomy,fin_repset):
    global repsetdic,taxdic,taxgendic,repsetIDlist,repgenlist,generaSeqIDdic
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
        genus = genus[-2]
        if genus[0:3] == "g__":
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
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if os.path.getsize(file) < 1:
            os.remove(file)
    return repsetdic,taxdic,taxgendic,repsetIDlist,repgenlist,generaSeqIDdic

def AlignToTree():
    #Muscle -> mini ITS alignments
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if file.endswith("_aligned.fasta"):
            continue
        if file.endswith("_seqs.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+muscledir+" -in "+inputname+" -out " +name+ "_aligned.fasta -quiet -maxiters 2 -diags1")
    #FastTree -> mini ITS trees
    for file in os.listdir(cwd):
        if file.endswith("_tree.nwk"):
            continue
        if file.endswith("_aligned.fasta"):
            inputname = str(file)
            file = file.split(".")
            name = file[0]
            os.system(""+ftdir+" -nt -quiet "+inputname+" > " +name+ "_tree.nwk")
 
    
def InsertITSin18S(backbone,hybridtree):
    #Open the 18S backbone tree file (newick format)
    with open (backbone, "r") as silvafile:
        finaltext=silvafile.read()
    #insert ITS-1 mini trees into 18S backbone
    cwd = os.getcwd()
    for file in os.listdir(cwd):
        if file.endswith("_tree.nwk"):
            genusname = str(file)
            genusname = genusname.split("_")
            genusname = genusname[2]
            str(genusname)
            with open (file, "r") as ITS1file:
                ITS1text=ITS1file.read().replace(";","")
            finaltext = finaltext.replace(genusname,ITS1text)
    fout = open(hybridtree,"w")
    finaltext = finaltext.replace(";","")
    finaltext = finaltext.replace("\n","")
    finaltext+= ";"
    fout.write(finaltext)
    fout.close()


#########################################
#  Directories & MANDATORY input files
#########################################
    
muscledir = "/Applications/muscle"                  #directory of MUSCLE
ftdir = "/Applications/./FastTree"                  #directory of FASTTREE
fin_taxonomy = "99_otu_taxonomy.txt"                #input file
fin_repset = "rep_set1.fna"                          #input file
silvadb = "SSURef_NR99_115_tax_silva_full_align_trunc.fasta"                             #input file

##  Name of ouput files (optional UNNECESSARY changes)
#########################################
fout = "FungiOnlySilvaAlignment.fasta"                       #name you want
fixedfastaname = "FungiOnlySilvaAlignment_fixed.fasta"      #name you want
backbone = "rep_phylo18Sbackbone.nwk"                        #name you want
hybridtree = "00_hybridtree.nwk"                             #name you want

##  call the functions
#########################################
tic()
MakeGeneraFastas(fin_taxonomy,fin_repset)
print "MakeGeneraFiles: "
toc()
ReduceSilva(silvadb,fout)
print "ReduceSilva: "
toc()
FixFastaFile(fout,fixedfastaname)
os.system(""+ftdir+" -nt -quiet "+fixedfastaname+" > "+backbone+"")
AlignToTree()
print "AlignToTree: "
toc()
InsertITSin18S(backbone,hybridtree)
print "InsertITSin18S: "
toc()

#These are indeed not even in the SILVA database
##Do not appear to be very abundant fungi (can calculate)
Warnings()

#os.system("filter_tree.py -i 082414_HybridTree.nwk -f rep_set1.fna -o 082414_HybridTree_pruned.nwk") 

######Uncomment these to delete extra files (leave commented for viewing genera files)
"""
for file in os.listdir(cwd):
    if file.endswith("_tree.nwk"):
        os.remove(file)
    if file.endswith("_aligned.fasta"):
        os.remove(file)
    if file.endswith("_seqs.fasta"):
        os.remove(file)
        
os.remove(backbone)
"""
