ghost-tree
==========

|Build Status| |Coverage Status|

ghost-tree is a bioinformatics tool that combines sequence data from two
genetic marker databases into one phylogenetic tree that can be used for
diversity analyses. One database is used as a backbone or scaffold because it
provides better phylogeny across all phyla, and the other database provides
finer taxonomic resolution.

For application to ITS, you don't need to install ghost-tree, but can use our
pre-built trees. The tree that you download needs to match the UNITE database
that you used (or plan to use) for your ITS analyses. The ``trees`` directory
contains pre-built reference phylogenetic trees for the `UNITE QIIME reference
files available here
<https://unite.ut.ee/repository.php>`_.

The most recent ghost-trees we've created are for the
`sh_qiime_release_s_30.12.2014 release of UNITE, available here
<https://unite.ut.ee/sh_files/sh_qiime_release_s_30.12.2014.zip>`_.
Depending on which of the fasta files you're using from that directory,
you'd use the corresponding ghost-tree listed below:

 * For ``sh_refs_qiime_ver6_97_30.12.2014.fasta``, use `ghosttree_UNITEv6_30.12.2014S_97_100clusters_052515.nwk <https://raw.githubusercontent.com/JTFouquier/ghost-tree/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_97_100clusters_052515.nwk>`_
 *  For ``sh_refs_qiime_ver6_97_30.12.2014.fasta``, use  `ghosttree_UNITEv6_30.12.2014S_97_80clusters_052515.nwk <https://github.com/JTFouquier/ghost-tree/raw/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_97_80clusters_052515.nwk>`_
 * For ``sh_refs_qiime_ver6_99_30.12.2014.fasta``, use `ghosttree_UNITEv6_30.12.2014S_99_100clusters_052515.nwk <https://raw.githubusercontent.com/JTFouquier/ghost-tree/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_99_100clusters_052515.nwk>`_
 *  For ``sh_refs_qiime_ver6_99_30.12.2014.fasta``, use  `ghosttree_UNITEv6_30.12.2014S_99_80clusters_052515.nwk <https://github.com/JTFouquier/ghost-tree/raw/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_99_80clusters_052515.nwk>`_
 * For ``sh_refs_qiime_ver6_dynamic_30.12.2014.fasta``, use `ghosttree_UNITEv6_30.12.2014S_dynamic_100clusters_052515.nwk <https://raw.githubusercontent.com/JTFouquier/ghost-tree/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_dynamic_100clusters_052515.nwk>`_
 *  For ``sh_refs_qiime_ver6_dynamic_30.12.2014.fasta``, use  `ghosttree_UNITEv6_30.12.2014S_dynamic_80clusters_052515.nwk <https://github.com/JTFouquier/ghost-tree/raw/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_dynamic_80clusters_052515.nwk>`_

Using ghost-tree.nwk files for your analyses:

To use the ghost-tree.nwk files in scripts such as
beta_diversity_through_plots.py in QIIME, you will need to filter your .biom
table so that it doesn't contain extra OTUs that will cause
beta_diversity_through_plots.py to fail. Note: We understand that QIIME isn't
the only downstream use for ghost-tree, but this has been a popular user
request.

This file `get_otus_from_ghost_tree.py <insertlink>`_ can be downloaded and
used to create a .txt file containing only the accession numbers from the
ghost-tree.nwk that you will be using for your diversity analyses.

You must have skbio installed to use `get_otus_from_ghost_tree.py`.
See: http://scikit-bio.org/

You will then use "ghost_tree_tips.txt" output file (containing the accession
numbers from the ghost-tree.nwk) to filter your .biom table so that it contains
only the OTUs that are in the ghost-tree.nwk that you are using.

The script, `filter_otus_from_otu_table.py <http://qiime.org/scripts/filter_otus_from_otu_table.html>`_
will filter your .biom table.

Use the required arguments in `filter_otus_from_otu_table.py` and also include
the following two arguments:
-e, --otu_ids_to_exclude_fp
(provide the text file containing OTU ids to exclude)
--negate_ids_to_exclude
(this will keep OTUs in otu_ids_to_exclude_fp, rather than discard them)

You should then have your filtered .biom table, a ghost-tree.nwk, and a mapping
file which will then allow you to use `beta_diversity_through_plots.py`
in QIIME.

**If** you are an experienced developer, or are interested in trying out the
ghost-tree tool via command line, then you will need to follow these
directions:

ghost-tree requires two external software tools to build a hybrid-tree or
the "ghost-tree":

MUSCLE (Version 3.8.31):
http://www.drive5.com/muscle/downloads.htm

FastTree (Version >2.1.7):
http://www.microbesonline.org/fasttree/#Install

To optionally (recommended) regroup the extension OTUs, ghost-tree requires
SUMACLUST:

SUMACLUST (Version 1.0.01):
https://git.metabarcoding.org/obitools/sumaclust/wikis/home/

Please install the software and make sure that it is in your PATH variable.
To test this you need to be able to type “muscle, fasttree or sumaclust” from
command line and see the “usage” or “help” documentation for each
software tool.

You will then need to "clone" the ghost-tree repository in order to download
all of the necessary files. Once you have cloned, you can then find the
`setup.py` and pip install it via "pip install -e ."

Typing the command `ghost-tree` will then display the ghost-tree help page
that provides command and subcommand help documentation.

This project is currently under active development. If you're interested in
contributing, please contact `@JTFouquier <https://github.com/JTFouquier>`__.

.. |Build Status| image:: https://travis-ci.org/JTFouquier/ghost-tree.svg?branch=master
   :target: https://travis-ci.org/JTFouquier/ghost-tree
.. |Coverage Status| image:: https://coveralls.io/repos/JTFouquier/ghost-tree/badge.png
   :target: https://coveralls.io/r/JTFouquier/ghost-tree
