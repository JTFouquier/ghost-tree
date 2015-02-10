import skbio
from skbio import BiologicalSequence


""" takes 97_OTUs (at "species" level) and clusters them at a lower level of
similarity using swarm.
"""


def preprocess_tip_sequences(species_level_otus_f):
    for seq in skbio.read(species_level_otus_f, format="fasta"):
        if "-" not in seq.id:
            new_seq_id = str(seq.id) + "_1"
        seq = BiologicalSequence(seq.sequence, id=new_seq_id)
        yield seq

# swarmdir = "/Users/jenniferfouquier/dev/ghost-tree/swarm"
