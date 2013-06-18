'''Check sequences for repeats that may impact cloning efficiency.'''

from collections import Counter
# TODO: implement version that works for circular DNA?


class FindRepeats(object):
    def __init__(self, dna_object, length):
        self.template = dna_object
        self.length = length

    def run(self):
        check = repeat_check(str(self.template), self.length)
        return check


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
    repeats = [(key, value) for key, value in counted.iteritems() if value > 1]
    return repeats
