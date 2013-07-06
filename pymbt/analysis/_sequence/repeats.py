'''
Check sequences for repeats that may impact cloning efficiency.

'''

from collections import Counter

# TODO: implement version that works for circular DNA


class Repeats(object):
    '''
    Count repeats in a DNA sequence.

    '''
    def __init__(self, dna, size):
        '''
        param dna_object: DNA sequence.
        type dna_object: DNA
        param size: size of repeat to detect (in bp).
        type size: int

        '''

        self.template = dna
        self.size = size

    def run(self):
        '''
        Execute function.

        '''
        check = repeats(self.template, self.size)
        return check


def repeats(sequence, size):
    '''
    Evaluate sequence repeats in a given sequence.

    :param seq: Input sequence.
    :type seq: str
    :param size: Size of the repeat to count.
    :type size: int

    '''
    sequence = str(sequence)

    n_mers = [sequence[i:i + size] for i in range(len(sequence) - size + 1)]
    counted = Counter(n_mers)
    repeats = [(key, value) for key, value in counted.iteritems() if value > 1]
    return repeats
