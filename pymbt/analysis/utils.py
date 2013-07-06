'''
Utils for analysis module.

'''

import pymbt.sequence


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
