"""
This file can be downloaded and used to create a .txt file containing only
the accession numbers from the ghost-tree.nwk that you plan to use for your
analyses.

You must have skbio installed. http://scikit-bio.org/
If you aren't familiar with skbio, make sure to check it out on its own, too!

You will then use "ghost_tree_tips.txt" output file containing the accession
numbers to filter your .biom table so that it contains only the OTUs that
are in the ghost-tree.nwk that you are using.

http://qiime.org/scripts/filter_otus_from_otu_table.html

Use the required arguments and the following two optional arguments:
-e, --otu_ids_to_exclude_fp
(provide the text file containing OTU ids to exclude)
--negate_ids_to_exclude
(this will keep OTUs in otu_ids_to_exclude_fp, rather than discard them)
"""
from skbio import TreeNode

ghosttree = TreeNode.read("ghost_tree_97_80clusters_from_alpha_release.nwk",
                          convert_underscores=False)  # your file goes here
output = open("ghost_tree_tips_underscore_fix.txt", "w")

for node in ghosttree.tips():
    output.write(str(node.name)+"\n")

output.close()
