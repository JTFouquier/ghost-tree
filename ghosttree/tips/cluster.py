import os


def preprocess_tip_sequences(tips_sequences_fh, similarity_threshold,
                             otu_formatted_fh):
    """Creates an OTU cluster formatted file from a .fasta sequence file.

    Creates an OTU cluster formatted file from a sequence file in .fasta
    format using Sumaclust (http://metabarcoding.org/sumatra). Sequence file
    can be from a database or from representative sequences specific to one's
    own dataset. OTU clusters will be created according to the threshold.
    The lower the threshold, the more "unidentified" sequences will be
    accounted for, but this will create larger "tips" groups that form less
    accurate multiple sequence alignments. Changing the threshold value to a
    number drastically different from the standard 97 or 99 OTU similarity
    will alter the resulting tree improving or reducing its utility in
    diversity analyses and should be performed with caution.

    Parameters
    __________
    tips_sequences_fh : filehandle
        Must be a .fasta formatted file from a database or representative
        sequences from a dataset. Usually these are 97 or 99 percent OTU
        sequences.

    similarity_threshold : float
        Value ranging from 0.00 to 1.00. Changing the threshold value to a
        number drastically different from the standard 97 or 99 OTU similarity
        will alter the resulting tree improving or reducing its utility in
        diversity analyses and should be performed with caution.

    otu_formatted_fh : filehandle
        Name of resulting OTU cluster formatted filehandle.

    """
    similarity_threshold = str(similarity_threshold)
    os.system("sumaclust -g -f -t "+similarity_threshold+" -O " +
              "" + otu_formatted_fh + " " + tips_sequences_fh + "")
