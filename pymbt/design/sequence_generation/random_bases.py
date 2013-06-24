'''
Generate random DNA or RNA sequences.

'''

import random
from pymbt import sequence


class RandomBases(object):
    '''
    Generate random DNA or RNA sequences.

    '''

    def __init__(self, size):
        '''
        :param size: Output sequence length.
        :type size: int

        '''

        self.size = size

    def __repr__(self):
        '''
        Printed representation of RandomBases object.

        '''
        return 'RandomBases generator for {} bases of DNA'.format(self.size)

    def run(self):
        '''
        Generate the sequence.

        '''

        random_seq = ''.join([random.choice('ATGC') for x in range(self.size)])
        return sequence.DNA(random_seq)
