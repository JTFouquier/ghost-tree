There have been a few reports in 2022 of incompatability with ghost-tree and newer releases of UNITE databases. As I have moved into a different role I am no longer actively maintaining this project. 

If anyone is interested in identifying the issues and creating a fix, you are welcome to. I will gladly give you a mention here. Thank you. 

ghost-tree
==========

ghost-tree is a bioinformatics tool that combines sequence data from two
genetic marker databases into one phylogenetic tree that can be used for
diversity analyses. One database is used as a "foundation tree" because it
provides better phylogeny across all phyla, and the other database provides
finer taxonomic resolution.

For application to ITS, you don't need to install *ghost-tree* but can use our
pre-built trees as the reference tree for your phylogenetic diversity analysis.

The tree that you download needs to match the UNITE database
that you used (or plan to use) for your ITS analyses. The most recent
ghost trees we've created are in our
`tree repository <https://github.com/JTFouquier/ghost-tree-trees>`_. This
repository contains pre-built reference phylogenetic trees for the
`UNITE QIIME reference database available here <https://unite.ut.ee/repository.php>`_.

.. image:: https://github.com/JTFouquier/q2-ghost-tree/blob/master/images/Picture1.png
Fig 1. Saliva (blue) and restroom (red) ITS sequences compared using A) binary
jaccard, B) unweighted UniFrac with Muscle aligned ITS sequences, and C)
unweighted UniFrac with a ghost tree created from ITS and 18S sequences.

Please cite our
`software publication in Microbiome <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0153-6>`_.

J. Fouquier, J.R. Rideout, E. Bolyen, J. Chase, A. Shiffer, D. McDonald, 
R. Knight, J.G. Caporaso, and S.T. Kelley. ghost-tree: creating hybrid-gene 
phylogenetic trees for diversity analysis. Microbiome. 
(February 2016) DOI: 10.1186/s40168-016-0153-6

Using ghost-tree.nwk files for your analyses:
=============================================

## QIIME 2

Directions are coming soon! You should be using Q2 for analysis now.
For now, please see `the q2-ghost-tree community tutorial
<https://github.com/JTFouquier/q2-ghost-tree/blob/master/QIIME2_community_tutorial.md>`_

## QIIME 1:
*please note that QIIME 1 is no longer officially supported. You should be
using QIIME 2 now*

To use the ghost-tree.nwk files in scripts, such as
`beta_diversity_through_plots.py
<http://qiime.org/scripts/beta_diversity_through_plots.html>`_
in QIIME, you will need to filter your .biom table so that it doesn't contain
extra OTUs that will cause `beta_diversity_through_plots.py
<http://qiime.org/scripts/beta_diversity_through_plots.html>`_ to fail.
Note: We understand that QIIME isn't the only downstream use for ghost-tree,
but this has been a popular user request.

This file, `get_otus_from_ghost_tree.py
<https://github.com/JTFouquier/ghost-tree/blob/master/helper_files/get_otus_from_ghost_tree.py>`_,
can be downloaded and used to create ``ghost_tree_tips.txt`` output file
containing only the accession numbers from the ghost-tree.nwk that you will
be using for your diversity analyses. You must have skbio installed to use
`get_otus_from_ghost_tree.py
<https://github.com/JTFouquier/ghost-tree/blob/master/helper_files/get_otus_from_ghost_tree.py>`_.
See `scikit-bio <http://scikit-bio.org/>`_ for install directions. scikit-bio
is very handy! You'll love it.

You will then use ``ghost_tree_tips.txt`` output file (containing the accession
numbers from the ghost-tree.nwk) to filter your .biom table so that it contains
only the OTUs that are in the ghost-tree.nwk that you are using.

The script, `filter_otus_from_otu_table.py
<http://qiime.org/scripts/filter_otus_from_otu_table.html>`_
will filter your .biom table.

Use the required arguments in `filter_otus_from_otu_table.py
<http://qiime.org/scripts/filter_otus_from_otu_table.html>`_ and also include
the following two arguments: `-e`, `--otu_ids_to_exclude_fp`
(provide the text file containing OTU ids to exclude) `--negate_ids_to_exclude`
(this will keep OTUs in otu_ids_to_exclude_fp, rather than discard them)

You should then have your filtered .biom table, a ghost-tree.nwk, and a mapping
file, which will then allow you to use `beta_diversity_through_plots.py
<http://qiime.org/scripts/beta_diversity_through_plots.html>`_
in QIIME.

Developers or to use ghost-tree via command line:
=================================================

If you are an experienced developer or are interested in trying out the
ghost-tree tool via command line, then you will need to follow the following
directions to start using *ghost-tree*.

ghost-tree requires two external software tools to build a hybrid-tree or
the "ghost-tree":

*Note: these three dependencies are all in Bioconda. You can install them using
`conda install the-software-tool -c bioconda`*

For additional information, see:

MUSCLE (Version 3.8.31):
http://www.drive5.com/muscle/downloads.htm

FastTree (Version >2.1.7):
http://www.microbesonline.org/fasttree/#Install

To optionally regroup the extension OTUs (recommended), ghost-tree requires
SUMACLUST:

SUMACLUST (Version 1.0.01):
https://git.metabarcoding.org/obitools/sumaclust/wikis/home/

Please install the software and make sure that it is in your PATH variable.
To test this you need to be able to type “muscle", "fasttree" or "sumaclust”
in the command line and see the corresponding “usage” or “help” documentation
for each software tool.

You will then need to "clone" the ghost-tree repository to download
all of the necessary files. You can then find the `setup.py` and install it via
"pip install -e ."

Typing the command `ghost-tree` will then display the *ghost-tree* help page
that provides command and subcommand help documentation.

You should also check out our "ipython notebook".  See the `ipynb install
directions here <http://ipython.org/install.html>`_.
For a detailed explanation of how to create your own ghost-tree.nwk
using the command line tool, see our `ghost-tree .ipynb workflow
<https://github.com/JTFouquier/ghost-tree/blob/master/workflow/ghost-tree_workflow.ipynb>`_.

This project is currently under active development and may evolve without
notice. So, please update your local repository, and check here for changes.
If you're interested in contributing, please contact
`@JTFouquier <https://github.com/JTFouquier>`_.

Help:
=====

If you have any trouble, please email jennietf@gmail.com. I am happy to help! :)

I am also interested in improving documentation, so please let me know if you
find errors or have suggestions!

.. |Build Status| image:: https://travis-ci.org/JTFouquier/ghost-tree.svg?branch=master
   :target: https://travis-ci.org/JTFouquier/ghost-tree
.. |Coverage Status| image:: https://coveralls.io/repos/JTFouquier/ghost-tree/badge.png
   :target: https://coveralls.io/r/JTFouquier/ghost-tree
