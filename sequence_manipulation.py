from string import maketrans
from string import translate

from pymbt.common_data import alphabets
from pymbt.common_data import codons 

# Simple reverse complement function
def reverse_complement(sequence):
    check_alphabet(sequence, material='dna')
    submat = maketrans('ATGCNatgcn', 'TACGNtacgn')
    sub_sequence = translate(sequence, submat)
    rev_sequence = sub_sequence[::-1]
    return rev_sequence


def check_alphabet(sequence, material='dna'):
    errs = {'dna': 'DNA', 'rna': 'RNA', 'pep': 'peptide'}

    '''Verifies that a given sequence is made only of DNA, RNA, or peptides'''
    if material is 'dna' or material is 'rna' or material is 'pep':
        alphabet = alphabets[material]
        err_msg = errs[material] 
    else:
        msg = 'Input material must be \'dna\', \'rna\', or \'pep\'.'
        raise ValueError(msg)
    for char in sequence:
        if char not in alphabet:
            raise ValueError('Sequence has a non-%s character' % err_msg)
    return sequence

def translate_seq(sequence):
    '''Input is DNA, output is peptide.
    This script expects only peptide sequence.
    If it doesn't find that, it errors'''
    # Split into 3-letter chunks
    sequence = check_alphabet(sequence, material='dna')
    # Make sure it's divisible by 3
    if (sum([int(x) for x in str(len(sequence))]) % 3) is not 0:
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
        aa = codons[codon]
        peptide.append(aa)
    peptide = ''.join(peptide)

    return peptide
