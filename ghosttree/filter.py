# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import skbio


def filter_positions(alignment_fh, maximum_gap_frequency,
                     maximum_position_entropy):
    """Filter gaps and high entropy positions from an alignment."""
    # (TODO)
    aln = skbio.alignment.TabularMSA(alignment_fh, 'phylip')

    aln = _filter_gap_positions(aln, maximum_gap_frequency)
    aln = _filter_high_entropy_positions(aln, maximum_position_entropy)
    return aln


def _filter_gap_positions(aln, maximum_gap_frequency):
    aln = aln.omit_gap_positions(maximum_gap_frequency)
    return aln


def _filter_high_entropy_positions(aln, maximum_position_entropy):
    entropies = aln.position_entropies(nan_on_non_standard_chars=False)
    positions_to_keep = []
    for position, entropy in enumerate(entropies):
        if entropy < maximum_position_entropy:
            positions_to_keep.append(position)
    aln = aln.subalignment(positions_to_keep=positions_to_keep)
    return aln
