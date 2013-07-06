'''
Helper functions for manipulating DNA, RNA, and peptide sequences.

'''

import re
from pymbt.data.common import ALPHABETS
from pymbt.data.common import COMPLEMENTS
# TODO: arguably this should not be imported - circular dependencies
import pymbt.sequence


def reverse_complement(sequence, material):
    '''
    Reverse complement a DNA sequence.

    :param sequence: Input sequence.
    :type sequence: str
    :param material: 'dna' or 'rna'.
    :type material: str

    '''

    complements = COMPLEMENTS[material]
    origin = complements[0]
    destination = complements[1]
    code = dict(zip(origin, destination))
    complemented = ''.join(code.get(k, k) for k in sequence)
    reverse_complemented = complemented[::-1]
    return reverse_complemented


def check_alphabet(sequence, material):
    '''
    Verify that a given string is made only of DNA, RNA, or peptide characters.

    :param sequence: DNA, RNA, or peptide sequence.
    :type sequence: str
    :param material: Input material - 'dna', 'rna', or 'pepide'.
    :type sequence: str

    '''

    errs = {'dna': 'DNA', 'rna': 'RNA', 'peptide': 'peptide'}
    if material == 'dna' or material == 'rna' or material == 'peptide':
        alphabet = ALPHABETS[material]
        err_msg = errs[material]
    else:
        msg = "Input material must be 'dna', 'rna', or 'peptide'."
        raise ValueError(msg)
    # This is a bottleneck for a lot of code.
    # First attempt with cython was slower. Could also try pypy.
    if re.search('[^' + alphabet + ']', sequence):
        raise ValueError('Sequence has a non-%s character' % err_msg)
    return sequence


def sequence_type(seq):
    '''
    Validates a DNA or RNA sequence instance.

    :param sequence_in: input DNA sequence.
    :type sequence_in: DNA
    :param material: 'dna' or 'rna'.
    :type material: str

    '''
    if isinstance(seq, pymbt.sequence.DNA):
        material = 'dna'
    elif isinstance(seq, pymbt.sequence.RNA):
        material = 'rna'
    elif isinstance(seq, pymbt.sequence.Peptide):
        material = 'peptide'
    else:
        raise Exception("Input was not a recognized pymbt.sequence object.")

    return material


def check_seq(seq, material):
    '''Do input checks / string processing.'''
    check_alphabet(seq, material)
    seq = seq.lower()
    return seq
