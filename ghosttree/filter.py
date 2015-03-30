
import skbio
from skbio.alignment import Alignment
from skbio.sequence import DNA

import numpy as np
# unit testing here is going to be difficult. Maybe move these into two
# functions
# I made three functions... and yes unit testing is a difficult here....
def filter_positions(alignment_fh, maximum_gap_frequency,
                     maximum_position_entropy):
    """Filter gaps and high entropy positions from an alignment."""
    aln = skbio.Alignment.read(alignment_fh)
    aln = _filter_gap_positions(aln, maximum_gap_frequency)
    aln = _filter_high_entropy_positions(aln, maximum_position_entropy)
    return aln

def _filter_gap_positions(aln, maximum_gap_frequency):
    aln = aln.omit_gap_positions(1.0 - np.finfo(float).eps)
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
