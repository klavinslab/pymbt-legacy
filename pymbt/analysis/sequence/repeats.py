'''Check sequences for repeats that may impact cloning efficiency.'''

from collections import Counter

from pymbt.sequence.utils import check_instance
# TODO: implement version that works for circular DNA?


class FindRepeats(object):
    '''
    Count repeats in a DNA sequence.

    '''
    def __init__(self, dna_object, size):
        '''
        param dna_object: DNA sequence.
        type dna_object: DNA
        param size: size of repeat to detect (in bp).
        type size: int

        '''

        self.template = dna_object
        check_instance(self.template)
        self.size = size

    def run(self):
        '''
        Execute function.

        '''
        check = _repeat_check(str(self.template), self.size)
        return check


def _repeat_check(sequence, size):
    '''
    Evaluate sequence repeats in a given sequence.

    :param seq: Input sequence.
    :type seq: str
    :param size: Size of the repeat to count.
    :type size: int

    '''

    n_mers = [sequence[i:i + size] for i in range(len(sequence) - size + 1)]
    counted = Counter(n_mers)
    repeats = [(key, value) for key, value in counted.iteritems() if value > 1]
    return repeats
