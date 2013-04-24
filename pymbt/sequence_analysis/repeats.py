'''Check sequences for repeats that may impact cloning efficiency.'''

from collections import Counter


def repeat_check(sequence, length):
    '''
    Evaluate sequence repeats in a given sequence.

    :param seq: Input sequence.
    :type seq: str
    :param length: Length of the repeat to count.
    :type length: int

    '''

    n_mers = [sequence[i:i + length] for i in range(len(sequence) - length)]
    counted = Counter(n_mers)
    return counted.most_common()
