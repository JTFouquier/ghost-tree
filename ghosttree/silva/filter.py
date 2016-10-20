# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import skbio
# (TODO) should I import the entire module or sections? We did both previously


def fungi_from_fasta(fasta_fh, accession_fh, taxonomy_fh):
    """Filter SILVA sequences to keep only fungi.

    Filters a fasta file of aligned or unaligned sequences to include only
    fungi. Only keeps sequences that have accession numbers that can be mapped
    to a fungal taxonomy string that ends at the genus rank.

    Parameters
    ----------
    fasta_fh : filehandle
        Fasta file of aligned or unaligned SILVA sequences. Each sequence
        identifier must be an accession number.
    accession_fh : filehandle
        A tab-separated file mapping accession numbers to a mapping number in
        `taxonomy_map`. This file should contain exactly two columns:
        accession number and mapping number.
    taxonomy_fh: filehandle
        A tab-separated file that identifes the taxonomy and rank of a mapping
        number in `accession_fh`. This file should contain exactly five
        columns beginning with taxonomy, mapping number and rank. The last two
        columns are ignored.

    Returns
    -------
    generator
        Yields ``skbio.BiologicalSequence`` objects.

    """
    accession_map = _parse_accession_map(accession_fh)
    taxonomy_map = _parse_taxonomy_map(taxonomy_fh)
    for seq in skbio.read(fasta_fh, format="fasta"):
        map_num = accession_map[seq.metadata['id']]
        if map_num in taxonomy_map:
            yield seq


def _parse_accession_map(accession_fh):
    accession_map = {}
    for line in accession_fh:
        accession, map_num = line.rstrip("\n").split("\t")
        if accession in accession_map:
            raise ValueError("Duplicate accession number %r" % accession)
        accession_map[accession] = map_num
    return accession_map


def _parse_taxonomy_map(taxonomy_fh):
    taxonomy_map = {}
    for line in taxonomy_fh:
        taxonomy, map_num, rank, _, _ = line.rstrip("\n").split("\t")
        if map_num in taxonomy_map:
            raise ValueError("Duplicate map number %r" % map_num)
        if rank == "genus" and "Fungi" in taxonomy:
            taxonomy_map[map_num] = taxonomy
    return taxonomy_map
