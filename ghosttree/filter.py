# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
from skbio import TabularMSA, DNA, RNA


def filter_positions(alignment_fh, maximum_gap_frequency,
                     maximum_position_entropy):
    """Filter gaps and high entropy positions from an alignment."""

    with alignment_fh:
        try:
            aln = TabularMSA.read(alignment_fh, constructor=DNA)
        except ValueError:
            alignment_fh.seek(0)
            aln = TabularMSA.read(alignment_fh, constructor=RNA)
    aln = _filter_gap_positions(aln, maximum_gap_frequency)
    aln = _filter_high_entropy_positions(aln, maximum_position_entropy)
    return aln


def _filter_gap_positions(aln, maximum_gap_frequency):

    aln_gap_frequencies = (aln.gap_frequencies(axis='sequence',
                                               relative=False) /
                           aln._seqs.count())
    aln_gap_frequencies_boolean = (aln_gap_frequencies <=
                                   maximum_gap_frequency)
    aln = aln.iloc[:, aln_gap_frequencies_boolean]

    return aln


def _filter_high_entropy_positions(aln, maximum_position_entropy):

    aln_entropies = aln.conservation(metric='inverse_shannon_uncertainty',
                                     degenerate_mode='nan', gap_mode='include')
    aln_entropies = 1 - aln_entropies
    aln = aln[:, (aln_entropies <= maximum_position_entropy)]

    return aln
