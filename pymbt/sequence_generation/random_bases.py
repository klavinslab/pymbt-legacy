'''Generate random DNA or RNA sequences.'''

import random


class RandomBases(object):
    '''Generate random DNA or RNA sequences.'''

    def __init__(self, size):
        '''
        :param size: Output sequence length.
        :type size: int

        '''

        self.size = size

    def __repr__(self):
        return 'RandomBases generator for {} bases of DNA'.format(self.size)

    def generate(self):
        '''Generate the sequence.'''

        random_seq = ''.join([random.choice('ATGC') for x in range(self.size)])
        return random_seq
