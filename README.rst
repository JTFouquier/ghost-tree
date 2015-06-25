ghost-tree
==========

|Build Status| |Coverage Status|

ghost-tree is a bioinformatics tool that combines sequence data from two
genetic marker databases into one phylogenetic tree that can be used for
diversity analyses. One database is used as a backbone or scaffold because it
provides better phylogeny across all phyla, and the other database provides
finer taxonomic resolution.

For application to ITS, you don't need to install ghost-tree, but can use our pre-built trees. The ``trees`` directory contains pre-built reference phylogenetic trees for the `UNITE QIIME reference files available here <https://unite.ut.ee/repository.php>`_.

The most recent ghost-trees we've created are for the `sh_qiime_release_s_30.12.2014 release of UNITE, available here <https://unite.ut.ee/sh_files/sh_qiime_release_s_30.12.2014.zip>`_. Depending on which of the fasta files you're using from that directory, you'd use the corresponding ghost-tree listed below:

 * For ``sh_refs_qiime_ver6_97_30.12.2014.fasta``, use `ghosttree_UNITEv6_30.12.2014S_97_100clusters_052515.nwk <https://raw.githubusercontent.com/JTFouquier/ghost-tree/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_97_100clusters_052515.nwk>`_
 *  For ``sh_refs_qiime_ver6_97_30.12.2014.fasta``, use  `ghosttree_UNITEv6_30.12.2014S_97_80clusters_052515.nwk <https://github.com/JTFouquier/ghost-tree/raw/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_97_80clusters_052515.nwk>`_
 * For ``sh_refs_qiime_ver6_99_30.12.2014.fasta``, use `ghosttree_UNITEv6_30.12.2014S_99_100clusters_052515.nwk <https://raw.githubusercontent.com/JTFouquier/ghost-tree/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_99_100clusters_052515.nwk>`_
 *  For ``sh_refs_qiime_ver6_99_30.12.2014.fasta``, use  `ghosttree_UNITEv6_30.12.2014S_99_80clusters_052515.nwk <https://github.com/JTFouquier/ghost-tree/raw/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_99_80clusters_052515.nwk>`_
 * For ``sh_refs_qiime_ver6_dynamic_30.12.2014.fasta``, use `ghosttree_UNITEv6_30.12.2014S_dynamic_100clusters_052515.nwk <https://raw.githubusercontent.com/JTFouquier/ghost-tree/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_dynamic_100clusters_052515.nwk>`_
 *  For ``sh_refs_qiime_ver6_dynamic_30.12.2014.fasta``, use  `ghosttree_UNITEv6_30.12.2014S_dynamic_80clusters_052515.nwk <https://github.com/JTFouquier/ghost-tree/raw/master/trees/ghost-trees_052515/ghosttree_UNITEv6_30.12.2014S_dynamic_80clusters_052515.nwk>`_




ghost-tree requires two external software tools to build a hybrid-tree or
the "ghost-tree":

MUSCLE (Version 3.8.31):
http://www.drive5.com/muscle/downloads.htm

FastTree (Version >2.1.7):
http://www.microbesonline.org/fasttree/#Install

To optionally (recommended) regroup the extension OTUs, ghost-tree requires
SUMACLUST:

SUMACLUST (Version 1.0.01):
http://metabarcoding.org/sumatra/wiki/download

Please install the software and make sure that it is in your PATH variable.
To test this you should be able to type “muscle, fasttree or sumaclust” from
command line and see the “usage” or “help” documentation for each
software tool.

This project is currently under active development. If you're interested in
contributing, please contact `@JTFouquier <https://github.com/JTFouquier>`__.

.. |Build Status| image:: https://travis-ci.org/JTFouquier/ghost-tree.svg?branch=master
   :target: https://travis-ci.org/JTFouquier/ghost-tree
.. |Coverage Status| image:: https://coveralls.io/repos/JTFouquier/ghost-tree/badge.png
   :target: https://coveralls.io/r/JTFouquier/ghost-tree
