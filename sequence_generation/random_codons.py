# TODO: implement 'dna' functionality

import random
from pymbt.common_data import codon_table
from pymbt.common_data import codon_freq_sc_nested
from pymbt.sequence_manipulation import check_alphabet


class RandomCodons:
    '''Provides a generator class for random DNA or RNA sequence.
         sequence: sequence for which to generate randomized codons.
         material: 'dna' for DNA or 'pep' for amino acid sequence.'''
    def __init__(self, sequence, freq_table = 'sc', threshold=0.0):
        self.sequence = check_alphabet(sequence, material='pep')
        # TODO: should generate this table from raw frequencies
        self.frequencies = codon_freq_sc_nested
        self.threshold = threshold

    def __repr__(self):
        return 'RandomCodons generator for %s' % self.sequence

    def generate(self):
        thresholded_codon_table = {}
        for key, value in self.frequencies.iteritems():
            average = 1.0 / len(value)
            new_value = [i for i, x in value.iteritems() if x/average >= self.threshold]
            if len(new_value) is 0:
                raise ValueError('The threshold has been set so high that it excludes all codons of a given amino acid.')
            thresholded_codon_table[key] = new_value
        coding_sequence = [random.choice(thresholded_codon_table[x]) for x in self.sequence]
        coding_sequence = ''.join(coding_sequence)
        return coding_sequence
