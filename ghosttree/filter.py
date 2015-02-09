
import skbio

# unit testing here is going to be difficult. Maybe move these into two
# functions
# make three functions ()
def filter_positions(alignment_fh, maximum_gap_frequency, maximum_position_entropy):
    """Filter gaps and high entropy positions from an alignment.

    """
    aln = skbio.Alignment.read(alignment_fh)
    aln = aln.omit_gap_positions(maximum_gap_frequency)
    aln = aln.omit_gap_positions(1.0 - np.finfo(float).eps)
    entropies = aln.position_entropies()
    positions_to_keep = []
    for position, entropy in enumerate(entropies):
        if entropy < maximum_position_entropy:
            positions_to_keep.append(position)
    aln = aln.subalignment(positions_to_keep=positions_to_keep)
    return aln
