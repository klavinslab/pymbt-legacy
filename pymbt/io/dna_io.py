# Get rid of this dependency?
from Bio import SeqIO

from pymbt import sequence


def read_dna(path, file_format):
    # This needs to be cleaned up - is only barely working for now
    seq = SeqIO.read(path, file_format)
    dna = sequence.DNA(seq.seq.tostring())
    dna.name = seq.name
    return dna
