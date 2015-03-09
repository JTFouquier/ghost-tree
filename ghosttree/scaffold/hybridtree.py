import re
import os

import skbio

from skbio import TreeNode
from skbio import Alignment
from skbio import read


def scaffold_tips_into_backbone(otu_table_fh, tips_taxonomy_fh, tips_seq_fh,
                                backbone_alignment_fh, ghost_tree_fp):
    """Inserts miniature trees into scaffold/backbone phylogenetic tree

    To Do:
    Need to create phylogenetic tree for the backbone tree.
    Need to align each genus' accession list and create phylogenetic tree.


    Parameters
    __________
    backbone_alignment_fh : filehandle
        Filepath containing aligned sequences from a genetic marker database
        in .fastq format.

    all_genus_dic

    Returns
    _______


    """
    global backbone_accession_genus_dic
    backbone_accession_genus_dic = {}
    all_genus_dic = _make_tips_genus_accession_dic(otu_table_fh,
                                                   tips_taxonomy_fh)
    skbio.write(_make_nr_backbone_alignment(backbone_alignment_fh,
                all_genus_dic),
                into="nr_backbone_alignment.fasta",
                format="fasta")
    backbone_tree = _make_backbone_tree("nr_backbone_alignment.fasta")
    all_tips_seqs = Alignment.read(tips_seq_fh)  # not alignment; needed
    # subset
    # for node in the backbone tree
    # Need to remove redundant seqs prior to making a tree
    for node in backbone_tree.tips():
        key_node, _ = str(node).split(":")
        key_node = backbone_accession_genus_dic[key_node]
        try:
            otu_seqs = all_tips_seqs.subalignment(all_genus_dic[key_node])
            mini_seq_file = open("mini_seq-gt.fasta", "w")
            mini_seq_file.write(str(otu_seqs))
            mini_seq_file.close()
            os.system(""+muscledir+" -in mini_seq-gt.fasta -out" +
                      " mini_alignment-gt.fasta -quiet -maxiters 2 -diags1")
            os.system(""+ftdir+" -nt -quiet mini_alignment-gt.fasta >" +
                      " mini_tree-gt.nwk")
            mini_tree = read("mini_tree-gt.nwk", format='newick',
                             into=TreeNode)
            node.append(mini_tree)  # insert mini trees
        except:
            continue
    ghost_tree_fp.write(str(backbone_tree))
    return str(backbone_tree).strip()


"""
complexity n^2; could be improved by having more genera than in tips file
currently contains only genera in tips file = n^2
Currently backbone is based on Accession number. Option to replace accession
with genus
"""


def _make_nr_backbone_alignment(backbone_alignment_fh, all_genus_dic):
    all_genus_list = all_genus_dic.keys()
    global backbone_accession_genus_dic
    backbone_accession_genus_dic = {}
    # currently n^2 This will take a while :(  need to go through genus dic
    # only need a backbone that contains seqs found in tips (genus_dic)
    for seq in skbio.read(backbone_alignment_fh, format="fasta"):
        try:
            for i in all_genus_list:
                if_case = (re.search(";" + i + ";", seq.description) or
                           re.search("g__" + i + ";", seq.description))
                if if_case:
                    all_genus_list.remove(i)
                    backbone_accession_genus_dic[seq.id] = i
                    yield seq
        except:
            pass


def _make_backbone_tree(in_name):
    os.system(""+ftdir+" -nt -quiet "+in_name+" > nr_backbone_tree.nwk")
    backbone_tree = TreeNode.read("nr_backbone_tree.nwk")
    return backbone_tree


def _make_tips_genus_accession_dic(otu_table_fh, tips_taxonomy_fh):
    """Find representative genus for each "tip cluster"

    This function takes in an OTU table file handle clustered at the user's
    desired percent similarity. This OTU table can be clustered once (i.e.)
    at 97 or 99 percent similarity) or can be from a secondary clustering step
    to capture more of the unidentified sequences so that they will be
    placed onto the final tree (see "ghost-tree group-tips"). Each OTU
    group will be XXXX

    Parameters
    __________
    otu_table_fh : filehandle
        The OTU table filehandle will be an OTU table of clusters where each
        line corresponds to a cluster of sequences represented by tab
        delimited accession numbers. OTUs were clustered based on a user's
        chosen percent similarity. First accession number can either be
        duplicated or not duplicated depending on OTU table format.
        example:
        A111\tA111\tA112
        A222\tA222\tA223
        A333\tA333\tA334
    tips_taxonomy_fh : filehandle
        The taxonomy file handle will be in a tab delimited file containing
        accession number and the corresponding taxonomy line. There are two
        columns, which are accession number and taxonomy line. The taxonomy
        line must be in the format =
        k__fungi;p__phylum;c__class;o__order;f__family;g__genus;s__species
        There are always two underscores following the taxonomy designation
        which is typical for "QIIME style" taxonomy lines (cite).
        example:
        A112\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Phoma;s__El
    modified_otu_table_fh : filehandle
        Table containing genus name and
        modified_otu_table_NR_fh : filehandle
        Needs to be non-redundant for genus, but contain groups of OTUs from
        different lines and clusters XXXXX

    Returns
    _______
    all_genus_dic : dict
    """
    accession_taxonomy_dic = _create_taxonomy_dic(tips_taxonomy_fh)
    all_genera_list = []
    global all_genus_dic
    all_genus_dic = {}
    for line in otu_table_fh:
        accession_list = line.strip().split("\t")  # line's accession list
        if accession_list[0] == accession_list[1]:
            del accession_list[0]  # remove the duplicate if there is one
        otu_genus_list = []  # changes for each OTU
        for i in accession_list:
            full_taxonomy_line = accession_taxonomy_dic[i]
            genus = full_taxonomy_line.split(";")
            genus = genus[-2]
            genus = genus[3:].capitalize()
            otu_genus_list.append(genus)
        most_common_genus = max(set(otu_genus_list), key=otu_genus_list.count)
        # genus winner for OTU in line
        all_genera_list.append(most_common_genus)  # genera for **ALL lines*
        if most_common_genus not in all_genus_dic:
            all_genus_dic[most_common_genus] = accession_list
        else:
            for i in accession_list:  # not efficient
                all_genus_dic[most_common_genus].append(i)
    otu_table_fh.close()
    return all_genus_dic


def _create_taxonomy_dic(tips_taxonomy_fh):
    accession_taxonomy_dic = {}
    for line in tips_taxonomy_fh:
        if "g__" not in line:
            raise ValueError("Taxonomy file must contain genera")
        accession, full_taxonomy_line = line.rstrip("\n").split("\t")
        accession = accession.strip()
        full_taxonomy_line = full_taxonomy_line.strip()
        accession_taxonomy_dic[accession] = full_taxonomy_line
    line = ""
    tips_taxonomy_fh.close()
    return accession_taxonomy_dic

muscledir = "/Applications/muscle"
ftdir = "/Applications/./FastTree"
