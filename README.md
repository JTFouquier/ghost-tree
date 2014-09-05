hybrid_tree
===========

This is a two genetic marker phylogenetic tree: hybrid tree

8/25/14: I combined many pieces of code -> 1 hybrid tree script 
  This complete script was tested many >10 times and a complete 18S/ITS tree that worked in UniFrac was produced
  
9/4/14: Of course there is a bug today. make_phylogeny.py works for 18S backbone, but FastTree won't make mini ITS trees

Bug:
Traceback (most recent call last):
  File "/macqiime/QIIME/bin/filter_tree.py", line 110, in <module>
    main()
  File "/macqiime/QIIME/bin/filter_tree.py", line 97, in main
    tree = DndParser(open(input_tree_fp,'U'))
  File "/macqiime/lib/python2.7/site-packages/cogent/parse/tree.py", line 156, in DndParser
    curr_node.Name = t
AttributeError: 'NoneType' object has no attribute 'Name'

9/5/14-- I think I know what I did... will fix soon






