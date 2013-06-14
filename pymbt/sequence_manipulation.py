'''Helper functions for manipulating DNA, RNA, and peptide sequences.'''

from string import maketrans
from pymbt.common_data import ALPHABETS
from pymbt.common_data import CODONS


def reverse_complement(sequence):
    '''
    Reverse complement a DNA sequence.

    :param sequence: DNA sequence.
    :type sequence: str

    '''

    check_alphabet(sequence, material='dna')
    submat = maketrans('ATGCNatgcn', 'TACGNtacgn')
    complemented = sequence.translate(submat)
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
    for char in sequence:
        if char not in alphabet:
            raise ValueError('Sequence has a non-{} character'.format(err_msg))
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
