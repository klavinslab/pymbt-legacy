'''Utils for analysis module.'''

from pymbt import sequence


def sequence_type(seq):
    '''Validates a pymbt.sequence data type.

    :param sequence_in: input DNA sequence.
    :type sequence_in: any

    '''
    if isinstance(seq, sequence.DNA):
        material = 'dna'
    elif isinstance(seq, sequence.RNA):
        material = 'rna'
    elif isinstance(seq, sequence.Peptide):
        material = 'peptide'
    else:
        raise Exception("Input was not a recognized pymbt.sequence object.")
    return material
