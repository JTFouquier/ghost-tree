import skbio
import re
from skbio import BiologicalSequence


def preprocess_tip_sequences(species_level_otus_fh):
    """Cluster OTUs in small groups for alignment and tree construction.

    Draft notes:

    Takes 97_OTUs (at "species" level) and clusters them at a lower level of
    similarity using swarm.

    Yields an accession number with a sequence that has been stripped of
    characters other than A,T,C or G

    _1 is a placeholder for the swarm package that requires "_" + abundance.
    But abundance of sequences can also be found in Biome tables.

    Arbitrarily replace non ATCG characters with "A" (this is up for debate,
    but only changes percent of As ~1 percent)

    Taxonomy file is a tab-delimited text file wih the first column being an
    accession number and the second column is the taxonomy ranking in the
    format

    k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Sebacinales;
    f__Sebacinaceae;g__unidentified;s__unculturedSebacina

    Parameters
    ----------
    species_level_otus_fh : filehandle
        Fasta file containing representative OTUs from a database or dataset.
        Each sequence identifier must be an accession number

    Returns
    -------
    generator
        Yields ``skbio.BiologicalSequence`` objects.

    """
    for seq in skbio.read(species_level_otus_fh, format="fasta"):
        seq = seq.upper()
        if "-" not in seq.id:
            new_seq_id = str(seq.id) + "_1"
        new_seq_sequence = re.sub('[^ATCG]', 'A', seq.sequence)
        seq = BiologicalSequence(new_seq_sequence, id=new_seq_id)
        yield seq


# muscledir = "/Applications/muscle"

# ftdir = "/Applications/./FastTree"
# swarmdir = "/Users/jenniferfouquier/dev/ghost-tree/swarm"

# os.system(""+muscledir+" -in "+inputname+" -out " +name+ "_aligned.fasta
# -quiet -maxiters 2 -diags1")

# os.system(""+ftdir+" -nt -quiet "+inputname+" > " +name+ "_tree.nwk")
