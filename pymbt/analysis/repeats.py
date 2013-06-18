'''Check sequences for repeats that may impact cloning efficiency.'''

from collections import Counter
# TODO: implement version that works for circular DNA?


class FindRepeats(object):
    def __init__(self, dna_object, size):
        self.template = dna_object
        self.size = size

    def run(self):
        check = repeat_check(str(self.template), self.size)
        return check


def repeat_check(sequence, size):
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
