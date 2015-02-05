import re
import os,sys

def ReduceSilva(fin,fout):
    fout1 = "intermediate.fasta"
    fin = open(fin,"U")
    fout1 = open(fout1,"w")
    fungi = 0
    for line in fin:
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
    fin.close()
    fout1.close()
    fin2 = open("intermediate.fasta","U")
    fout2 = "intermediate2.fasta"
    fout2 = open(fout2,"w")
    genuslist = []
    alignmentseq = ""
    genusname = ""
    for line in fin2:
        if re.search("-------",line):
            alignmentseq = line
        if re.search(">",line):
            line = re.sub("_Eukaryota","",line)
            line = line.split(";")
            genusorfamily = line[-2]
            genusorfamily = genusorfamily[-4:]
            if genusorfamily == "ceae":
                genusname = line[-1]
                genusname = re.split("_",genusname)
                genusname = genusname[0]
            if genusorfamily != "ceae":
                genusname = line[-2]
            genusname = ">" + genusname
            if genusname not in genuslist:
                fout2.write(alignmentseq)
                fout2.write("\n")
                fout2.write(genusname)
                fout2.write("\n")
                genuslist.append(genusname)
    fout2.write(alignmentseq)               
    fin2.close()
    fout2name = fout2.name
    fout2.close()
    os.remove("intermediate.fasta")
    os.system("fix_arb_fasta.py -f "+fout2name+" -o "+fout+"")
    os.remove("intermediate2.fasta")

fin = "SSURef_NR99_115_tax_silva_full_align_trunc.fasta"
fout = "FungiOnlySilvaAlignment.fasta"

###STEP 1
ReduceSilva(fin,fout)

os.system("fix_arb_fasta.py -f FungiOnlySilvaAlignment.fasta -o FungiOnlySilvaAlignment_fixed.fasta")
os.system("make_phylogeny.py -i FungiOnlySilvaAlignment_fixed.fasta -o rep_phylo18Sbackbone.nwk") 

###STEP 2

#Way bulky -> need to simplify with dictionaries:
#accession #, taxonomy & fasta sequence

f1in = "99_otu_taxonomy.txt"

f1in = open(f1in,"U")
glist = []
for line in f1in:
    line = line[8:]
    line = line.strip()
    line = line.split(";")
    g = line[5]
    if g not in glist:
        glist.append(g)
f1in.close() 

def SelectTaxonomy(f1in,taxon):
    f1in = open(f1in,"U")
    f1out = open(taxon+"_smalltaxonomy.txt","w")
    for line in f1in:
        if re.search(taxon,line):
            f1out.write(line)
    f1in.close()
    f1out.close()

#Example from TAXONOMY IDENTIFICATION FILE! Ex: Pick out all of the cladosporium seqs
#GQ458030	k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Capnodiales;f__Mycosphaerellaceae;g__Cladosporium;s__Cladosporium_cladosporioides

def IDlistfromRepSet(f1in,repset,taxon):
    f1in = open(f1in,"U")
    taxonNumList = []
    for line in f1in:
        line = line[0:8]
        taxonNumList.append(line)
    f1in.close()
    f2in = open(repset,"U")
    global seqIDsfromRepSetList
    seqIDsfromRepSetList = []
    for line in f2in:
        idinrepset = line[1:9]
        for k in taxonNumList:
            if k == idinrepset:
                seqIDsfromRepSetList.append(line)
    f2in.close()

def ExtractMySeqs(myrepset,myseqIDs,taxon):
    myseqIDlist = []
    myseqFastaFile = open(taxon+"_seqs.fasta","w") 
    myseqIDs = seqIDsfromRepSetList
    myrepset = open(myrepset,"U") 
    count=0
    nicdic={}
    nicval=""
    for line in myrepset:
                if re.match(">",line):
                        if len(nicval)>0:
                                nicdic[nickey]=nicval
                                count+=1
                        seqname=line
                        seqname=seqname.strip()
                        nickey=seqname
                        nicval=""
                else:
                        nicval+=line.strip()
    myrepset.close()
    for i in range(4):
        value = nicdic.popitem()
    for line in myseqIDs:
        line = line[0:]
        line = line.strip()
        myseqIDlist.append(line)
    IDlist = []
    error = 0
    for i in myseqIDlist:
        try: k = reverse_lookup(nicdic,nicdic[i])
        except:
            #print "warning ",error
            error+=1
            continue
        if k not in IDlist:
            out="" + k + "\n" + nicdic[i] + "\n"
            myseqFastaFile.write(out)
            IDlist.append(k)
        else:
            continue
    myseqFastaFile.close()

def reverse_lookup(nicdic, v):
    for k in nicdic:
        if nicdic[k] == v:
            return k

seqIDsfromRepSetList = []
for i in glist:
    taxon = i
    myrepset = "rep_set1.fna"
    SelectTaxonomy("99_otu_taxonomy.txt", taxon)
    IDlistfromRepSet(taxon+"_smalltaxonomy.txt",myrepset,taxon)
    myseqIDs = seqIDsfromRepSetList
    if myseqIDs == []:
        os.remove(taxon+"_smalltaxonomy.txt")
        continue
    ExtractMySeqs(myrepset,myseqIDs,taxon)
    os.remove(taxon+"_smalltaxonomy.txt")    

cwd = os.getcwd()
for file in os.listdir(cwd):
    if os.path.getsize(file) < 1:
        os.remove(file)

#Muscle -> mini ITS alignments
cwd = os.getcwd()
for file in os.listdir(cwd):
    if file.endswith("_aligned.fasta"):
        continue
    if file.endswith("_seqs.fasta"):
        inputname = str(file)
        file = file.split(".")
        name = file[0]
        os.system("muscle -in "+inputname+" -out " +name+ "_aligned.fasta -quiet -maxiters 2 -diags1")

#FastTree -> mini ITS trees
cwd = os.getcwd()
for file in os.listdir(cwd):
    if file.endswith("_tree.nwk"):
        continue
    if file.endswith("_aligned.fasta"):
        inputname = str(file)
        file = file.split(".")
        name = file[0]
        os.system("FastTree -nt -quiet "+inputname+" > " +name+ "_tree.nwk")
        
#Open the 18S backbone tree file (newick format)
fin= "rep_phylo18Sbackbone.nwk"
with open (fin, "r") as silvafile:
    finaltext=silvafile.read()

#insert ITS-1 mini trees into 18S backbone
for file in os.listdir(cwd):
    if file.endswith("_tree.nwk"):
        genusname = str(file)
        genusname = genusname.split("_")
        genusname = genusname[2]
        str(genusname)
        with open (file, "r") as ITS1file:
            ITS1text=ITS1file.read().replace(";","")
        finaltext = finaltext.replace(genusname,ITS1text)

fout = open("082414_HybridTree.nwk","w")
finaltext = finaltext.replace(";","")
finaltext = finaltext.replace("\n","")
finaltext+= ";"
fout.write(finaltext)
fout.close()

os.system("filter_tree.py -i 082414_HybridTree.nwk -f rep_set1.fna -o 082414_HybridTree_pruned.nwk") 

#Uncomment these to delete extra files (leave commented for viewing genera files)
"""
for file in os.listdir(cwd):
    if file.endswith("_tree.nwk"):
        os.remove(file)
    if file.endswith("_aligned.fasta"):
        os.remove(file)
    if file.endswith("_seqs.fasta"):
        os.remove(file)
        
os.remove("rep_phylo18Sbackbone.nwk")
"""
