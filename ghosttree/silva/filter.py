import skbio


def fungi_from_fasta(fasta_fh, accession_fh, taxonomy_fh):
    accession_map = _parse_accession_map(accession_fh)
    taxonomy_map = _parse_taxonomy_map(taxonomy_fh)
    for seq in skbio.read(fasta_fh, format="fasta"):
        map_num = accession_map[seq.id]
        if map_num in taxonomy_map:
            yield seq


def _parse_accession_map(accession_fh):
    accession_map = {}
    for line in accession_fh:
        accession, map_num = line.strip().split("\t")
        if accession in accession_map:
            raise ValueError("Duplicate accession number %r" % accession)
        accession_map[accession] = map_num
    return accession_map


def _parse_taxonomy_map(taxonomy_fh):
    taxonomy_map = {}
    for line in taxonomy_fh:
        split_line = line.strip().split("\t")
        if len(split_line) == 3:
            taxonomy, map_num, rank = split_line
            if map_num in taxonomy_map:
                raise ValueError("Duplicate map number %r" % map_num)
            if rank == "genus" and "Fungi" in taxonomy:
                taxonomy_map[map_num] = taxonomy
    return taxonomy_map
