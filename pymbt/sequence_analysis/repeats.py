from collections import Counter


def repeat_check(seq, n):
    '''
    Evaluate sequence repeats in a given sequence.

    :param seq: Input sequence.
    :type seq: str.
    :param n: Repeat size to find.
    :type n: int.

    '''
    n_mers = [seq[i:i + n] for i in range(len(seq) - n)]
    counted = Counter(n_mers)
    return counted.most_common()
