# Get rid of this dependency?
import os
from Bio import SeqIO
from pymbt import sequence


def read_dna(path, file_format):
    # This needs to be cleaned up - is only barely working for now
    seq = SeqIO.read(path, file_format)
    dna = sequence.DNA(seq.seq.tostring())
    dna.name = seq.name
    return dna


def read_dnas(dirpath, file_format):
    '''
    Read .seq results files from a dir.

    :param dirpath: Path to directory containing sequencing files.
    :type dirpath: str
    :param file_format: Format of file. e.g. genbank, fasta, abi, seq
    :type file_format: str

    '''

    seq_paths = [x for x in os.listdir(dirpath) if x.endswith('.seq')]
    abi_paths = [x for x in os.listdir(dirpath) if x.endswith('.abi')]
    abi_paths += [x for x in os.listdir(dirpath) if x.endswith('.ab1')]
    seq_seqs = [read_dna(dirpath + x, 'fasta') for x in seq_paths]
    abi_seqs = [read_dna(dirpath + x, 'abi') for x in abi_paths]
    sequences = seq_seqs + abi_seqs
    return sequences
