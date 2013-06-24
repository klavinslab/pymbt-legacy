'''Helper functions for manipulating DNA, RNA, and peptide sequences.'''

import math
import re
from pymbt.data.common import ALPHABETS
from pymbt.data.common import COMPLEMENTS
from pymbt.data.common import CODONS
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


def check_alphabet(sequence, material='dna'):
    '''
    Verify that a given string is made only of DNA, RNA, or peptide characters.

    :param sequence: DNA, RNA, or peptide sequence.
    :type sequence: str
    :param material: Input material - 'dna', 'rna', or 'pep'.
    :type sequence: str

    '''

    errs = {'dna': 'DNA', 'rna': 'RNA', 'pep': 'peptide'}
    if material == 'dna' or material == 'rna' or material == 'pep':
        alphabet = ALPHABETS[material]
        err_msg = errs[material]
    else:
        msg = 'Input material must be \'dna\', \'rna\', or \'pep\'.'
        raise ValueError(msg)
    # This is a bottleneck for a lot of code.
    # First attempt with cython was slower. Could also try pypy.
    if re.search('[^' + alphabet + ']', sequence):
        raise ValueError('Sequence has a non-%s character' % err_msg)
    return sequence


def translate_seq(sequence):
    '''
    Translate a DNA sequence into peptide sequence.

    :param sequence: DNA sequence.
    :type sequence: str

    '''

    # Split into 3-letter chunks
    sequence = check_alphabet(sequence, material='dna')
    # Make sure it's divisible by 3
    if (sum([int(x) for x in str(len(sequence))]) % 3) != 0:
        raise Exception('Input was not a valid coding sequence')
    n_codons = len(sequence) / 3
    # Make a list of codons (3 char strings)
    sequence = list(sequence)
    peptide = []
    for i in range(n_codons):
        base_1 = sequence.pop(0)
        base_2 = sequence.pop(0)
        base_3 = sequence.pop(0)
        codon = ''.join(base_1 + base_2 + base_3).upper()
        amino_acid = CODONS[codon]
        peptide.append(amino_acid)
    while peptide[-1] == '*':
        peptide.pop()
    peptide = ''.join(peptide)

    return peptide


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
    elif isinstance(seq, pymbt.sequence.DNA):
        material = 'rna'
    else:
        raise Exception("Input was not a DNA or RNA object.")

    return material


def check_seq(seq, material):
    '''Do input checks / string processing.'''
    check_alphabet(seq, material=material)
    seq = seq.lower()
    return seq


def check_inv(pattern):
    '''
    Check whether pattern is palindrome.
    :param pattern: pattern to test.
    :type pattern: str

    '''
    p_len = len(pattern)
    wing = int(math.floor(p_len / 2))
    if p_len % 2 != 0:
        l_wing = pattern[0:wing + 1]
        r_wing = pattern[wing:]
    else:
        l_wing = pattern[0: wing]
        r_wing = pattern[wing:]
    if l_wing == reverse_complement(r_wing, 'dna'):
        return True
    else:
        return False
