ghost-tree
==========

|Build Status| |Coverage Status|

ghost-tree is a bioinformatics tool that combines sequence data from two
genetic marker databases into one phylogenetic tree that can be used for
diversity analyses. One database is used as a backbone or scaffold because it
provides better phylogeny across all phyla, and the other database provides
finer taxonomic resolution.

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
