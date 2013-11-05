'''Resection (need a new name!) - exonuclease activity.'''


def five_resect(dna, n_bases):
    """Remove bases from 5' end of top strand.

    :param dna: Sequence to resect.
    :type dna: pymbt.sequence.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int
    :returns: DNA sequence resected at the 5' end by n_bases.
    :rtype: pymbt.sequence.DNA

    """
    new_instance = dna.copy()
    new_top = '-' * min(len(dna), n_bases) + str(dna)[n_bases:]
    new_instance._sequence = new_top
    new_instance._remove_end_gaps()
    if n_bases >= len(dna):
        new_instance._sequence = ''.join(['-' for i in range(len(dna))])
        new_instance.stranded = 'ss'
    return new_instance


def three_resect(dna, n_bases):
    '''Remove bases from 3' end of top strand.

    :param dna: Sequence to resect.
    :type dna: pymbt.sequence.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int
    :returns: DNA sequence resected at the 3' end by n_bases.
    :rtype: pymbt.sequence.DNA

    '''
    new_instance = dna.copy()

    new_top = str(dna)[:-n_bases] + '-' * min(len(dna), n_bases)
    new_instance._sequence = new_top
    new_instance._remove_end_gaps()
    if n_bases >= len(dna):
        new_instance._sequence = ''.join(['-' for i in range(len(dna))])
        new_instance.stranded = 'ss'
    return new_instance
