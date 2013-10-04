'''Utils for analysis module.'''
from pymbt import sequence


def sequence_type(seq):
    '''Validates a pymbt.sequence data type.

    :param sequence_in: input DNA sequence.
    :type sequence_in: any
    :returns: The material - 'dna', 'rna', or 'peptide'.
    :rtype: str
    :raises: ValueError

    '''
    if isinstance(seq, sequence.DNA):
        material = 'dna'
    elif isinstance(seq, sequence.RNA):
        material = 'rna'
    elif isinstance(seq, sequence.Peptide):
        material = 'peptide'
    else:
        raise ValueError("Input was not a recognized pymbt.sequence object.")
    return material
