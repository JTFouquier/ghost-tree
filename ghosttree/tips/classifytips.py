# from collections import Counter


# this step is all after doing Sumaclust
def find_rep_otu_genus(otu_table_fh, taxonomy_fh, modified_otu_table_fh,
                       modified_otu_table_NR_fh):
    accession_taxonomy_dic = _create_taxonomy_dic(taxonomy_fh)
    all_genera_list = []
    otu_genus_dic = {}
    for line in otu_table_fh:
        accession_list = line.strip().split("\t")
        del accession_list[0]
        otu_genus_list = []  # changes for each OTU
        for i in accession_list:
            full_taxonomy_line = accession_taxonomy_dic[i]
            genus = full_taxonomy_line.split(";")
            genus = genus[-2]
            genus = genus[3:]
            otu_genus_list.append(genus)
        # genus_count_dic = Counter(otu_genus_list)
        most_common_genus = max(set(otu_genus_list), key=otu_genus_list.count)
        all_genera_list.append(most_common_genus)
        accession_list_str = "\t".join(accession_list)
        modified_otu_table_fh.write(genus+"\t"+accession_list_str+"\n")
        if most_common_genus not in otu_genus_dic:
            otu_genus_dic[most_common_genus] = accession_list
        else:
            for i in accession_list:  # not efficient
                otu_genus_dic[most_common_genus].append(i)
    tuple_dic = otu_genus_dic.items()
    for i in tuple_dic:
        i = list(i)
        modified_otu_table_NR_fh.write(str(i)+"\n")
    otu_table_fh.close()
    modified_otu_table_fh.close()
    modified_otu_table_NR_fh.close()
    print otu_genus_dic
    return otu_genus_dic


def _create_taxonomy_dic(taxonomy_fh):
    line = ""
    accession_taxonomy_dic = {}
    for line in taxonomy_fh:
        accession, full_taxonomy_line = line.rstrip("\n").split("\t")
        accession = accession.strip()
        full_taxonomy_line = full_taxonomy_line.strip()
        accession_taxonomy_dic[accession] = full_taxonomy_line
    line = ""
    taxonomy_fh.close()
    return accession_taxonomy_dic

"""
needs more arguments... testing purposes only
find_rep_otu_genus("miniotus.txt", "minitaxonomy.txt",
                   "modifiedotu021615.txt")
find_rep_otu_genus("fullITSOTUmap97_clustered80.txt", "97_otu_taxonomy.txt",
                   "modifiedotularge021615.txt")
"""
