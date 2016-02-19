ghost-tree
==========

|Build Status| |Coverage Status|

ghost-tree is a bioinformatics tool that combines sequence data from two
genetic marker databases into one phylogenetic tree that can be used for
diversity analyses. One database is used as a "foundation tree" because it
provides better phylogeny across all phyla, and the other database provides
finer taxonomic resolution.

For application to ITS, you don't need to install ghost-tree, but can use our
pre-built trees. The tree that you download needs to match the UNITE database
that you used (or plan to use) for your ITS analyses. The ``trees`` directory
contains pre-built reference phylogenetic trees for the `UNITE QIIME reference
files available here
<https://unite.ut.ee/repository.php>`_.

The most recent ghost-trees we've created are in our tree repository.

Using ghost-tree.nwk files for your analyses:
=============================================

To use the ghost-tree.nwk files in scripts, such as
`beta_diversity_through_plots.py <http://qiime.org/scripts/beta_diversity_through_plots.html>`_
in QIIME, you will need to filter your .biom table so that it doesn't contain
extra OTUs that will cause `beta_diversity_through_plots.py <http://qiime.org/scripts/beta_diversity_through_plots.html>`_ to fail.
Note: We understand that QIIME isn't the only downstream use for ghost-tree,
but this has been a popular user request.

This file, `get_otus_from_ghost_tree.py <https://github.com/JTFouquier/ghost-tree/blob/master/helper_files/get_otus_from_ghost_tree.py>`_,
can be downloaded and used to create ``ghost_tree_tips.txt`` output file
containing only the accession numbers from the ghost-tree.nwk that you will
be using for your diversity analyses. You must have skbio installed to use `get_otus_from_ghost_tree.py <https://github.com/JTFouquier/ghost-tree/blob/master/helper_files/get_otus_from_ghost_tree.py>`_.
See `scikit-bio <http://scikit-bio.org/>`_ for install directions. scikit-bio
is very handy! You'll love it.

You will then use ``ghost_tree_tips.txt`` output file (containing the accession
numbers from the ghost-tree.nwk) to filter your .biom table so that it contains
only the OTUs that are in the ghost-tree.nwk that you are using.

The script, `filter_otus_from_otu_table.py <http://qiime.org/scripts/filter_otus_from_otu_table.html>`_
will filter your .biom table.

Use the required arguments in `filter_otus_from_otu_table.py <http://qiime.org/scripts/filter_otus_from_otu_table.html>`_ and also include
the following two arguments:
-e, --otu_ids_to_exclude_fp
(provide the text file containing OTU ids to exclude)
--negate_ids_to_exclude
(this will keep OTUs in otu_ids_to_exclude_fp, rather than discard them)

You should then have your filtered .biom table, a ghost-tree.nwk, and a mapping
file, which will then allow you to use `beta_diversity_through_plots.py <http://qiime.org/scripts/beta_diversity_through_plots.html>`_
in QIIME.

If you had any trouble, please email jennietf@gmail.com.
I am also interested in improving documentation, so please let me know if you
find errors or have suggestions! Thank you!

Developers:
===========

If you are an experienced developer, or are interested in trying out the
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

You should also really check out our "ipython notebook."  See `ipynb install directions here <http://ipython.org/install.html>`_.
For a detailed explanation on how to create your own ghost-tree.nwk
using the command line tool, see our `ghost-tree .ipynb workflow <https://github.com/JTFouquier/ghost-tree/blob/master/workflow/ghost-tree_workflow.ipynb>`_.

This project is currently under active development and may evolve so please
update your local repository and check here for changes. If you're interested in
contributing, please contact `@JTFouquier <https://github.com/JTFouquier>`__.

.. |Build Status| image:: https://travis-ci.org/JTFouquier/ghost-tree.svg?branch=master
   :target: https://travis-ci.org/JTFouquier/ghost-tree
.. |Coverage Status| image:: https://coveralls.io/repos/JTFouquier/ghost-tree/badge.png
   :target: https://coveralls.io/r/JTFouquier/ghost-tree
