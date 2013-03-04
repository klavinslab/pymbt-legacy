import random


class RandomBases:
    '''Provides a generator class for random DNA or RNA sequence'''
    def __init__(self, n):
        self.n = n

    def __repr__(self):
        return 'RandomBases generator for %s bases of DNA') % (self.n

    def generate(self):
        random_seq = ''.join([random.choice('ATGC') for x in range(self.n)])
        return random_seq
