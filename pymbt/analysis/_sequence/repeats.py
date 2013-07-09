'''
Check sequences for repeats that may impact cloning efficiency.

'''

from collections import Counter

# TODO: implement version that works for circular DNA


def repeats(seq, size):
    '''
    Evaluate sequence repeats in a given sequence.

    :param seq: Input sequence.
    :type seq: pymbt.sequence.DNA or pymbt.sequence.RNA
    :param size: Size of the repeat to count.
    :type size: int

    '''
    seq = str(seq)

    n_mers = [seq[i:i + size] for i in range(len(seq) - size + 1)]
    counted = Counter(n_mers)
    # No one cares about patterns that appear once, so exclude them
    repeats = [(key, value) for key, value in counted.iteritems() if value > 1]

    return repeats
