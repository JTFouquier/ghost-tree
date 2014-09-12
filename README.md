hybrid_tree
===========
<p>
This is a two genetic marker phylogenetic tree: hybrid tree

Installation instructions:<p>

Need:<p>
  1)repset of ITS seqs (I used UNITE)<p>
  2)taxonomy file<p>
  --above are in Google Drive (will share)<p>
  3)Silva 18S file:<p>
  http://www.arb-silva.de/fileadmin/silva_databases/release_115/Exports/SSURef_NR99_115_tax_silva_full_align_trunc.fasta.tgz
  (very large, full alignment, non-redundant SSU, has ALL eukaryotes which is unnecessary)
  Current release is 119 and 30GB!? why big change (JF look into)   
<p>
For now, this can be run in command line in QIIME in the working directory:  "python make_hybrid_tree.py"
<p>
It uses Muscle and Fasttree from QIIME (and other Qiime scripts)
<p>
And that is it.
<p>
<p>
<p>
Known issues/questions for later:
  9/7/14<p>
    JF --> use dictionary for organizing OTUs into genera files (IN PROGRESS via wip (work in progress))<p>
    JF --> make root to tip calculation/warnings. <p> 
    JF --> remove dependency on Qiime... fix fasta files, make_phylogeny, etc. (IN PROGRESS) <p>
    Q: Do people want the option to view their files (the seq file, alignment or tree). i.e. they want to look at the OTUs in Phoma or Cladosporium.  (can currently uncomment or comment this feature) <p>
    Q: SILVA is very messy. It is not organized as well as UNITE. i.e. sometimes when looking for a genus (by slicing the 2nd to last item in the list) it will grab a family.  So I added a check for that.  <p>
    Q: Warnings?  When your most abundant ITS genera are not in the Silva tree? That could mess up the tree, correct? <p> 










