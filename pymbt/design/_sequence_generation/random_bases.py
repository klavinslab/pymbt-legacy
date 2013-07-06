'''
Generate a random DNA sequence.

'''

import random
from pymbt import sequence


class RandomDNA(object):
    '''
    Generate a random DNA sequence.

    '''

    def __init__(self, length):
        '''
        :param length: Output sequence length.
        :type length: int

        '''

        self.length = length

    def __repr__(self):
        '''
        Printed representation of RandomBases object.

        '''
        return 'RandomBases generator for {} bases of DNA'.format(self.size)

    def choose_base(self):
        return sequence.DNA(random.choice('ATGC'))

    def run(self):
        '''
        Generate the sequence.

        '''

        random_seq = ''.join([self.choose_base() for i in range(self.length)])
        return sequence.DNA(random_seq)
