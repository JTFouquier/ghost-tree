# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import re
import os
import subprocess
import sys
import tempfile

import skbio
import pandas as pd


def extensions_onto_foundation(otu_file_fh, extension_taxonomy_fh,
                               extension_seq_fh,
                               foundation_fh,
                               ghost_tree_fp,
                               graft_level, foundation_taxonomy):
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

    foundation_fh : filehandle
        File containing EITHER pre-aligned sequences from a genetic marker
        database in .fasta format OR a newick tree. This file refers to the
        "foundation" of the ghost-tree.

        .fasta contains accession numbers *and* taxonomy labels.

        .nwk tree is a tree with accession numbers. MUST supply a foundation
        taxonomy if using a tree as a foundation. Pass this via the
        --foundation-taxonomy option. see --help for details

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

    graft_level, graft_letter = _graft_functions(graft_level)

    process = subprocess.Popen("muscle", shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    std_output, std_error = process.communicate()
    if re.search("command not found", str(std_error)):
        print("muscle, multiple sequence aligner, is not found. "
              "Is it installed? Is it in your path?")

    process = subprocess.Popen("fasttree", shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    std_output, std_error = process.communicate()
    if re.search("command not found", str(std_error)):
        print("fasttree, phylogenetic tree builder, is not found. "
              "Is it installed? Is it in your path?")
    os.mkdir(ghost_tree_fp)
    extension_genus_accession_list_dic = \
        _extension_genus_accession_dict(otu_file_fh, extension_taxonomy_fh,
                                        graft_level)

    sniffer_results = skbio.io.sniff(foundation_fh)[0]

    if sniffer_results == 'newick':
        if foundation_taxonomy is None:
            sys.exit("ghost-tree error: You must provide a foundation "
                     "taxonomy if using a foundation tree. Pass the taxonomy "
                     "file using the '--foundation-taxonomy' flag.")
        foundation_tree = \
            _make_nr_foundation_newick(foundation_fh,
                                       extension_genus_accession_list_dic,
                                       graft_letter, foundation_taxonomy)

    if sniffer_results == 'fasta':
        nr_foundation_alignment = \
            _make_nr_foundation_alignment(foundation_fh,
                                          extension_genus_accession_list_dic,
                                          graft_letter)
        skbio.io.write(nr_foundation_alignment,
                       into=ghost_tree_fp +
                       "/nr_foundation_alignment_gt.fasta",
                       format="fasta")
        foundation_tree, all_std_error = \
            _make_foundation_tree(ghost_tree_fp +
                                  "/nr_foundation_alignment_gt.fasta",
                                  str(std_error), ghost_tree_fp)

    seqs = list(skbio.io.read(extension_seq_fh, format='fasta'))
    all_std_error = ''

    with tempfile.TemporaryDirectory() as tmp:

        for node in foundation_tree.tips():
            key_node = node.name
            key_node = foundation_accession_genus_dic[key_node]
            _make_mini_otu_files(key_node, extension_genus_accession_list_dic,
                                 seqs, tmp)
            process = subprocess.Popen("muscle -in " + tmp +
                                       "/mini_seq_gt.fasta -out " + tmp +
                                       "/mini_alignment_gt.fasta -quiet" +
                                       " -maxiters 2 -diags1", shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            std_output, std_error = process.communicate()
            process = subprocess.Popen("fasttree -nt -quiet " +
                                       tmp + "/mini_alignment_gt.fasta >" +
                                       tmp + "/mini_tree_gt.nwk", shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            std_output, std_error = process.communicate()
            all_std_error += "FastTree warnings for genus " + key_node + \
                             " are:\n" + str(std_error) + "\n"
            mini_tree = skbio.io.read(tmp + "/mini_tree_gt.nwk",
                                      format='newick', into=skbio.TreeNode)
            node.extend(mini_tree.root_at_midpoint().children[:])
    # print('GRAFT LEVEL: ', graft_letter)
    # print(foundation_tree.ascii_art())
    ghost_tree_nwk = open(ghost_tree_fp + "/ghost_tree.nwk", "w")
    ghost_tree_nwk.write(str(foundation_tree))
    ghost_tree_nwk.close()
    _make_accession_id_file(ghost_tree_fp)
    return str(foundation_tree).strip(), all_std_error


def _make_accession_id_file(ghost_tree_fp):
    ghosttree = skbio.io.read(ghost_tree_fp + "/ghost_tree.nwk",
                              format='newick', into=skbio.TreeNode)
    output = open(ghost_tree_fp + "/ghost_tree_extension_accession_ids.txt",
                  "w")
    for node in ghosttree.tips():
        output.write(str(node.name) + "\n")
    output.close()


def _make_mini_otu_files(key_node, extension_genus_accession_list_dic, seqs,
                         tmp):
    keep = extension_genus_accession_list_dic[key_node]

    output_file = open(tmp + "/mini_seq_gt.fasta", "w")

    for seq in seqs:
        if seq.metadata['id'] in keep:
            fasta_format = ">" + seq.metadata['id'] + "\n" + str(seq) + "\n"
            output_file.write(fasta_format)
    output_file.close()


def _extension_genus_accession_dict(otu_file_fh, extension_taxonomy_fh,
                                    graft_level):
    """Find representative genus for each "extension tree cluster" """
    accession_taxonomy_dic = _create_taxonomy_dict(extension_taxonomy_fh,
                                                   graft_level)
    all_genera_in_extension_list = []
    extension_genus_accession_list_dic = {}
    for line in otu_file_fh:
        accession_list = line.strip().split("\t")
        if accession_list[0] == accession_list[1]:
            del accession_list[0]
        otu_genus_list = []
        for i in accession_list:
            taxonomy = accession_taxonomy_dic[i]
            otu_genus_list.append(taxonomy)
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


def _create_taxonomy_dict(extension_taxonomy_fh, graft_level):

    accession_id_list = []
    taxonomy_list = []
    for line in extension_taxonomy_fh:
        # Qiime2 specific
        if re.search("Feature ID", line):
            continue
        splitline = line.split('\t')
        accession = splitline[0].strip()
        taxonomy_line = splitline[1].strip()
        accession_id_list.append(accession)
        taxonomy_list.append(taxonomy_line)

    ds_taxonomy = pd.Series(taxonomy_list)
    accession_taxonomy_dict = _collapse_taxa_line(accession_id_list,
                                                  ds_taxonomy, graft_level)
    return accession_taxonomy_dict


def _make_nr_foundation_newick(foundation_fh,
                               extension_genus_accession_list_dic,
                               graft_letter, foundation_taxonomy):
    global foundation_accession_genus_dic
    foundation_accession_genus_dic = {}
    all_genus_list = list(extension_genus_accession_list_dic.keys())
    foundation_tree = skbio.io.read(foundation_fh, format='newick',
                                    into=skbio.TreeNode,
                                    convert_underscores=False)

    foundation_unique_accessions = []
    for line in foundation_taxonomy:
        splitline = line.split('\t')
        accession = splitline[0].strip()
        foundation_taxonomy = splitline[1].strip()
        for graft_taxa in all_genus_list:
            if_case = (re.search(";" + graft_taxa.lower() + ";",
                                 foundation_taxonomy.lower()) or
                       re.search(graft_letter + "__" + graft_taxa.lower() +
                                 ";", foundation_taxonomy.lower()) or
                       re.search(";" + graft_taxa.lower(),
                                 foundation_taxonomy.lower()))

            if if_case:
                all_genus_list.remove(graft_taxa)
                foundation_accession_genus_dic[accession] = graft_taxa
                foundation_unique_accessions.append(accession)

    sheared_tree = foundation_tree.shear(foundation_unique_accessions)
    return sheared_tree


def _make_nr_foundation_alignment(foundation_fh,
                                  extension_genus_accession_list_dic,
                                  graft_letter):
    all_genus_list = list(extension_genus_accession_list_dic.keys())
    global foundation_accession_genus_dic
    foundation_accession_genus_dic = {}
    for seq in skbio.io.read(foundation_fh, format="fasta"):
        for i in all_genus_list:

            if_case = (re.search(";" + i + ";", seq.metadata['description']) or
                       re.search(graft_letter + "__" + i + ";",
                                 seq.metadata['description']) or
                       re.search(";" + i, seq.metadata['description']))
            if if_case:
                all_genus_list.remove(i)
                foundation_accession_genus_dic[seq.metadata['id']] = i
                seq.metadata['description'] = i
                yield seq


def _make_foundation_tree(in_name, all_std_error, ghost_tree_fp):
    process = subprocess.Popen("fasttree -nt -quiet " + in_name + "" +
                               " > " +
                               ghost_tree_fp +
                               "/nr_foundation_tree_gt.nwk",
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    std_output, std_error = process.communicate()
    all_std_error += "Error log for ghost-tree:\n\n\nSome genera may not " \
                     "contain any errors, so the genus is listed as a " \
                     "placeholder\n\n"
    all_std_error += "FastTree warnings for the foundation_tree are:\n" + \
                     str(std_error) + "\n"

    foundation_tree = skbio.io.read(ghost_tree_fp +
                                    "/nr_foundation_tree_gt.nwk",
                                    format='newick', into=skbio.TreeNode)
    foundation_tree.root_at_midpoint()
    return foundation_tree, all_std_error


def _collapse_taxa_line(accession_ids: list, taxonomy: pd.Series,
                        graft_level: int) -> pd.DataFrame:
    """This function was copied and modified from q2-taxa from Qiime2."""
    # Assemble the taxonomy data
    max_observed_level = _get_max_level(taxonomy)
    if graft_level > max_observed_level:
        raise ValueError('Requested level of %d is larger than the maximum '
                         'level available in taxonomy data (%d).' %
                         (graft_level, max_observed_level))

    accession_taxonomy_dict = {}
    index_count = 0
    for tax in taxonomy:
        new_taxa = _collapse(tax, graft_level)
        accession_taxonomy_dict[accession_ids[index_count]] = new_taxa
        index_count += 1

    return accession_taxonomy_dict


def _get_max_level(taxonomy):
    """This function was copied from q2-taxa from Qiime2."""
    return taxonomy.apply(lambda x: len(x.split(';'))).max()


def _collapse(tax, level):
    """This function was copied and modified from q2-taxa from Qiime2."""
    tax = [x.strip() for x in tax.split(';')]
    taxa = tax[:level][-1].split('__')[1].capitalize()
    return taxa


def _graft_functions(graft_level):
    graft_letter = graft_level
    graft_level_map = {'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6}
    graft_level = graft_level_map[graft_letter]
    return graft_level, graft_letter
