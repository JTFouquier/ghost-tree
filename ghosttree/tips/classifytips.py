# from collections import Counter


# Needs to address OTU table errors (duplicates)
# this step is all after doing Sumaclust
def find_rep_otu_genus(otu_table_fh, taxonomy_fh, modified_otu_table_fh,
                       modified_otu_table_NR_fh):
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
    taxonomy_fh : filehandle
        The taxonomy file handle will be in a tab delimited file containing the
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
    accession_taxonomy_dic = _create_taxonomy_dic(taxonomy_fh)
    all_genera_list = []
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
        accession_list_str = "\t".join(accession_list)
        modified_otu_table_fh.write(genus+"\t"+accession_list_str+"\n")
        if most_common_genus not in all_genus_dic:
            all_genus_dic[most_common_genus] = accession_list
        else:
            for i in accession_list:  # not efficient
                all_genus_dic[most_common_genus].append(i)
    tuple_dic = all_genus_dic.items()
    for i in tuple_dic:
        i = list(i)
        modified_otu_table_NR_fh.write(str(i)+"\n")
    otu_table_fh.close()
    modified_otu_table_fh.close()
    modified_otu_table_NR_fh.close()
    print all_genus_dic
    return all_genus_dic


def _create_taxonomy_dic(taxonomy_fh):
    line = ""
    accession_taxonomy_dic = {}
    for line in taxonomy_fh:
        if "g__" not in line:
            raise ValueError("Taxonomy file must contain genera %r" %
                             line)
        accession, full_taxonomy_line = line.rstrip("\n").split("\t")
        accession = accession.strip()
        full_taxonomy_line = full_taxonomy_line.strip()
        accession_taxonomy_dic[accession] = full_taxonomy_line
    line = ""
    taxonomy_fh.close()
    return accession_taxonomy_dic
