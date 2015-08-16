# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import re
import os
import shutil
import subprocess

import skbio

from skbio import TreeNode
from skbio.alignment import SequenceCollection
from skbio import read


def extensions_onto_foundation(otu_file_fh, extension_taxonomy_fh,
                               extension_seq_fh,
                               foundation_alignment_fh,
                               ghost_tree_fp):
    """Combines two genetic databases into one phylogenetic tree.

    Some genetic databases provide finer taxonomic resolution,
    but high sequence variability causes poor multiple sequence alignments
    (these are the "extension trees"). Other databases provide high quality
    phylogenetic information (hence it is used as the "foundation"), but poor
    taxonomic resolution. This script combines two genetic databases into
    one phylogenetic tree in .nwk format, taking advantage of the benefits
    of both databases, but allowing sequencing to be performed using the
    "extension trees" primer set.

    Parameters
    __________
    otu_file_fh : filehandle
        Tab-delimited text file containing OTU clusters in rows containing
        accession numbers only. Format can be 1) where the accession number
        is in the first column with only one column or 2) it can contain
        accession numbers clustered in tab-delimited rows containing more
        accession numbers, which are part of that OTU cluster (as in output of
        "ghost-tree group-extensions"). This file refers to the "extension
        trees". File references to sequence reads or sample numbers/names are
        not valid here. This is not an OTU .biom table.

    extension_taxonomy_fh : filehandle
        Tab-delimited text file related to "extension trees" wih the 1st
        column being an
        accession number (same accession numbers in otu_file_fh and
        extension_taxonomy_fh) and the 2nd column is the taxonomy ranking in
        the following format:
        k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Sebacinales;
        f__Sebacinaceae;g__unidentified;s__Sebacina

    extension_seq_fh : filehandle
        The .fasta formated sequences for the "extension trees" genetic
        dataset. Sequence identifiers are the accession numbers. These
        accession numbers are the same as in the otu_file_fh and
        extension_taxonomy_fh.

    foundation_alignment_fh : filehandle
        File containing pre-aligned sequences from a genetic marker database
        in .fasta format. This file refers to the "foundation" of the
        ghost-tree. Contains accession numbers and taxonomy labels.

    ghost_tree_fp : folder
        Output folder contains files including:
        a) The Newick formatted ghost-tree, which is the final output of the
           ghost-tree tool. This is a phylogenetic tree designed for
           downstream diversity analyses.
        b) Accession IDs from the ghost-tree.nwk file that you can use for
           downstream analyses tools
        c) log error file (this is an optional file that you can have if you
           type '--stderr')
    """
    global foundation_accession_genus_dic  # needs global assignment for flake8
    foundation_accession_genus_dic = {}
    std_output, std_error = "", ""
    process = subprocess.Popen("muscle", shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    std_output, std_error = process.communicate()
    if re.search("command not found", std_error):
        print "muscle, multiple sequence aligner, is not found. Is it" \
              " installed? Is it in your path?"
    process = subprocess.Popen("fasttree", shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    std_output, std_error = process.communicate()
    std_output, std_error = "", ""
    if re.search("command not found", std_error):
        print "fasttree, phylogenetic tree builder, is not found. Is it" \
              " installed? Is it in your path?"
    os.mkdir("tmp")
    os.mkdir(ghost_tree_fp)
    extension_genus_accession_list_dic = \
        _extension_genus_accession_dic(otu_file_fh,
                                       extension_taxonomy_fh)
    skbio.write(_make_nr_foundation_alignment(foundation_alignment_fh,
                extension_genus_accession_list_dic),
                into=ghost_tree_fp + "/nr_foundation_alignment_gt.fasta",
                format="fasta")
    foundation_tree, all_std_error = _make_foundation_tree(ghost_tree_fp + "/nr_foundation_alignment_gt.fasta",
                                                           std_error, ghost_tree_fp)
    seqs = SequenceCollection.read(extension_seq_fh)
    for node in foundation_tree.tips():
        key_node, _ = str(node).split(":")
        key_node = foundation_accession_genus_dic[key_node]
        try:
            _make_mini_otu_files(key_node, extension_genus_accession_list_dic,
                                 seqs)
            process = subprocess.Popen("muscle -in tmp/mini_seq_gt.fasta" +
                                       " -out" +
                                       " tmp/mini_alignment_gt.fasta -quiet" +
                                       " -maxiters 2 -diags1", shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            std_output, std_error = process.communicate()
            process = subprocess.Popen("fasttree -nt -quiet" +
                                       " tmp/mini_alignment_gt.fasta >" +
                                       " tmp/mini_tree_gt.nwk", shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            std_output, std_error = process.communicate()
            all_std_error += "FastTree warnings for genus "+key_node+" are:\n" + std_error + "\n"
            mini_tree = read("tmp/mini_tree_gt.nwk", format='newick',
                             into=TreeNode)
            node.extend(mini_tree.root_at_midpoint().children[:])
        except:
            continue
    shutil.rmtree("tmp")
    ghost_tree_nwk = open(ghost_tree_fp + "/ghost_tree.nwk", "w")
    ghost_tree_nwk.write(str(foundation_tree))
    ghost_tree_nwk.close()
    _make_accession_id_file(ghost_tree_fp)
    return str(foundation_tree).strip(), all_std_error


def _make_accession_id_file(ghost_tree_fp):
    ghosttree = TreeNode.read(ghost_tree_fp + "/ghost_tree.nwk")
    output = open(ghost_tree_fp + "/ghost_tree_extension_accession_ids.txt", "w")
    for node in ghosttree.tips():
        output.write(str(node.name)+"\n")
    output.close()


def _make_mini_otu_files(key_node, extension_genus_accession_list_dic, seqs):
    keep = extension_genus_accession_list_dic[key_node]
    output_file = open("tmp/mini_seq_gt.fasta", "w")
    for seq in seqs:
        if seq.id in keep:
            fasta_format = ">"+seq.id+"\n"+str(seq)+"\n"
            output_file.write(fasta_format)
    output_file.close()
    return fasta_format


def _extension_genus_accession_dic(otu_file_fh, extension_taxonomy_fh):
    """Find representative genus for each "extension tree cluster" """
    accession_taxonomy_dic = _create_taxonomy_dic(extension_taxonomy_fh)
    all_genera_in_extension_list = []
    extension_genus_accession_list_dic = {}
    for line in otu_file_fh:
        accession_list = line.strip().split("\t")
        if accession_list[0] == accession_list[1]:
            del accession_list[0]
        otu_genus_list = []
        for i in accession_list:
            full_taxonomy_line = accession_taxonomy_dic[i]
            genus = full_taxonomy_line.split(";")
            genus = genus[-2]
            genus = genus[3:].capitalize()
            otu_genus_list.append(genus)
        most_common_genus = max(set(otu_genus_list), key=otu_genus_list.count)
        if most_common_genus == "Unidentified":
            otu_genus_set = set(otu_genus_list)
            otu_genus_set.remove("Unidentified")
            try:
                most_common_genus = max(otu_genus_set,
                                        key=otu_genus_list.count)
            except:
                pass
        all_genera_in_extension_list.append(most_common_genus)
        if most_common_genus not in extension_genus_accession_list_dic:
            extension_genus_accession_list_dic[most_common_genus] = \
                accession_list
        else:
            for i in accession_list:
                extension_genus_accession_list_dic[most_common_genus].append(i)
    otu_file_fh.close()
    return extension_genus_accession_list_dic


def _create_taxonomy_dic(extension_taxonomy_fh):
    accession_taxonomy_dic = {}
    for line in extension_taxonomy_fh:
        if "g__" not in line:
            raise ValueError("Taxonomy file must contain genera")
        accession, full_taxonomy_line = line.rstrip("\n").split("\t")
        accession = accession.strip()
        full_taxonomy_line = full_taxonomy_line.strip()
        accession_taxonomy_dic[accession] = full_taxonomy_line
    line = ""
    extension_taxonomy_fh.close()
    return accession_taxonomy_dic


def _make_nr_foundation_alignment(foundation_alignment_fh,
                                  extension_genus_accession_list_dic):
    all_genus_list = extension_genus_accession_list_dic.keys()
    global foundation_accession_genus_dic
    foundation_accession_genus_dic = {}
    for seq in skbio.read(foundation_alignment_fh, format="fasta"):
        try:
            for i in all_genus_list:
                if_case = (re.search(";" + i + ";", seq.description) or
                           re.search("g__" + i + ";", seq.description))
                if if_case:
                    all_genus_list.remove(i)
                    foundation_accession_genus_dic[seq.id] = i
                    yield seq
        except:
            pass


def _make_foundation_tree(in_name, all_std_error, ghost_tree_fp):
    process = subprocess.Popen("fasttree -nt -quiet "+in_name+"" +
                               " > "+ghost_tree_fp+"/nr_foundation_tree_gt.nwk", shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    std_output, std_error = process.communicate()
    all_std_error += "Error log for ghost-tree:\n\n\nSome genera may not contain " + \
                     "any errors, so the genus is listed as a placeholder\n\n"
    all_std_error += "FastTree warnings for the foundation_tree are:\n" + std_error + "\n"
    foundation_tree = TreeNode.read(ghost_tree_fp + "/nr_foundation_tree_gt.nwk")
    return foundation_tree, all_std_error
