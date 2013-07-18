'''Resection (need a new name!) - exonuclease activity.'''


def five_resect(dna, n_bases):
    '''Remove bases from 5' end of top strand.

    :param dna: Sequence to resect.
    :type dna: pymbt.sequence.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int

    '''
    new_instance = dna.copy()
    new_top = '-' * min(len(dna.top), n_bases) + dna.top[n_bases:]
    new_instance.top = new_top
    new_instance.remove_end_gaps()
    if n_bases >= len(dna):
        new_instance = dna.set_stranded('ss')
    return new_instance


def three_resect(dna, n_bases):
    '''Remove bases from 3' end of top strand.

    :param dna: Sequence to resect.
    :type dna: pymbt.sequence.DNA
    :param n_bases: Number of bases cut back.
    :type n_bases: int

    '''
    new_instance = dna.copy()

    new_top = dna.top[:-n_bases] + '-' * min(len(dna.top), n_bases)
    new_instance.top = new_top
    new_instance.remove_end_gaps()
    if n_bases >= len(dna):
        new_instance = dna.set_stranded('ss')
    return new_instance
