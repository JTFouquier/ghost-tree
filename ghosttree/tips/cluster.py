import skbio
import re
from skbio import BiologicalSequence


""" takes 97_OTUs (at "species" level) and clusters them at a lower level of
similarity using swarm.

Yields an accession number with a sequence that has been stripped of any
resricted characters

_1 is a placeholder for the swarm package that requires "_" + abundance.
But abundance of sequences can also be found in OTU tables.

Arbitrarily replace non ATCG characters with "A" (this is questionable)
"""


def preprocess_tip_sequences(species_level_otus_f):
    for seq in skbio.read(species_level_otus_f, format="fasta"):
        seq = seq.upper()
        if "-" not in seq.id:
            new_seq_id = str(seq.id) + "_1"
        new_seq_sequence = re.sub('[^ATCG]', 'A', seq.sequence)
        seq = BiologicalSequence(new_seq_sequence, id=new_seq_id)
        yield seq

# swarmdir = "/Users/jenniferfouquier/dev/ghost-tree/swarm"
