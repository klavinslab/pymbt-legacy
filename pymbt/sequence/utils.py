'''Helper functions for manipulating DNA, RNA, and peptide sequences.'''
import re
from pymbt.data.common import ALPHABETS
from pymbt.data.common import COMPLEMENTS


def reverse_complement(seq, material):
    '''Reverse complement a DNA or RNA sequence.

    :param seq: Input sequence.
    :type seq: str
    :param material: 'dna' or 'rna'.
    :type material: str

    '''
    complements = COMPLEMENTS[material]
    origin = complements[0]
    destination = complements[1]
    code = dict(zip(origin, destination))
    complemented = ''.join(code.get(base) for base in seq)
    reverse_complemented = complemented[::-1]
    return reverse_complemented


def check_alphabet(seq, material):
    '''Verify that a given string is valid DNA, RNA, or peptide characters.

    :param seq: DNA, RNA, or peptide sequence.
    :type seq: str
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
    # This is a bottleneck when modifying sequence - hence the run_checks
    # optional parameter in sequence objects..
    # First attempt with cython was slower. Could also try pypy.
    if re.search('[^' + alphabet + ']', seq):
        raise ValueError('Sequence has a non-%s character' % err_msg)


def process_seq(seq, material):
    '''Validate and process sequence inputs.

    :param seq: input sequence
    :type seq: pymbt.sequence.{DNA, RNA, Peptide}
    :param material: DNA, RNA, or peptide
    :type: str
    '''
    check_alphabet(seq, material)
    seq = seq.lower()
    return seq


def palindrome(seq):
    '''Test whether a sequence is palindrome.

    :param seq: Sequence to analyze (DNA or RNA).
    :type seq: pymbt.sequence.DNA or pymbt.sequence.RNA

    '''
    seq_len = len(seq)
    if seq_len % 2 == 0:
        # Sequence has even number of bases, can test non-overlapping seqs
        wing = seq_len / 2
        l_wing = seq[0: wing]
        r_wing = seq[wing:]
        if l_wing == r_wing.reverse_complement():
            return True
        else:
            return False
    else:
        # Sequence has odd number of bases and cannot be a palindrome
        return False
