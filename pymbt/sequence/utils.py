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


# TODO: Split up dna/rna conversion from translation?
def convert_sequence(seq, from_material, to_material):
    '''
    Translate a DNA sequence into peptide sequence.

    :param seq: DNA or RNA sequence.
    :type seq: DNA or RNA
    :param from_material: material to convert ('rna', or 'dna')
    :type from_material: str
    :param to_material: material to which to convert ('rna', 'dna', or
                        'peptide').
    :type to_material: str

    '''

    if from_material == 'dna' and to_material == 'rna':
        # Can't transcribe a gap
        if pymbt.sequence.DNA('-') in seq:
            raise ValueError('Cannot transcribe gapped DNA')
        # Convert DNA chars to RNA chars
        origin = ALPHABETS['dna'][:-1]
        destination = ALPHABETS['rna']
        code = dict(zip(origin, destination))
        converted = ''.join(code.get(str(k), str(k)) for k in seq)
        # Instantiate RNA object
        converted = pymbt.sequence.RNA(converted)
    elif from_material == 'rna' and to_material == 'dna':
        # Convert RNA chars to DNA chars
        origin = ALPHABETS['rna']
        destination = ALPHABETS['dna'][:-1]
        code = dict(zip(origin, destination))
        converted = ''.join(code.get(str(k), str(k)) for k in seq)
        # Instantiate DNA object
        converted = pymbt.sequence.DNA(converted)
    elif from_material == 'rna' and to_material == 'peptide':
        # Make a list for easier processing
        seq_list = list(str(seq))

        # Convert to peptide until stop codon is found.
        converted = []
        while True:
            if len(seq_list) >= 3:
                base_1 = seq_list.pop(0)
                base_2 = seq_list.pop(0)
                base_3 = seq_list.pop(0)
                codon = ''.join(base_1 + base_2 + base_3).upper()
                amino_acid = CODONS[codon]
                # Stop when stop codon is found
                if amino_acid == '*':
                    break
                converted.append(amino_acid)
            else:
                break
        converted = ''.join(converted)
        converted = pymbt.sequence.Peptide(converted)
    else:
        msg1 = 'Conversion from '
        msg2 = '{0} to {1} is not supported.'.format(from_material,
                                                     to_material)
        raise ValueError(msg1 + msg2)

    return converted


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
