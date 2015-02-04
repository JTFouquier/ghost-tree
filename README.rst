hybrid_tree
===========
<p>
This is a two genetic marker phylogenetic tree: hybrid tree

Installation instructions:<p>  

For BOTH scripts, you need 3 input files:<p>
  1)repset of ITS seqs<p>
  2)taxonomy file<p>
  --above are in Google Drive (will share)<p>
  3)Silva 18S file:<p>
  http://www.arb-silva.de/fileadmin/silva_databases/release_115/Exports/SSURef_NR99_115_tax_silva_full_align_trunc.fasta.tgz
  (very large, full alignment, non-redundant SSU, has ALL eukaryotes which is unnecessary)
<p>

There are two scripts, one relies on Qiime (make_hybrid_tree.py) and the other one (skbio_hybrid_tree_wip.py) uses skbio with FastTree and Muscle installed locally

"python make_hybrid_tree.py" can be run in command line in QIIME in the working directory where your 3 files are; It uses Muscle and Fasttree from within QIIME (and other Qiime scripts)
<p>
<p>
skbio_hybrid_tree_wip.py is the transition away from Qiime dependency; CHANGE input files at the bottom of python script & change directory of FastTree and Muscle (you must download and install locally)
<p>
<p>
