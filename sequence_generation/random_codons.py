import random
from pymbt.common_data import codon_table
from pymbt.common_data import codon_freq_sc_nested
from pymbt.sequence_manipulation import check_alphabet
from pymbt.sequence_manipulation import translate_seq


class RandomCodons:
    '''Provides a generator class for random DNA or RNA sequence.
         sequence: sequence for which to generate randomized codons.
         material: 'dna' for DNA or 'pep' for amino acid sequence.'''
    def __init__(self,
                 sequence,
                 material='pep',
                 freq_table='sc',
                 threshold=0.0):
        self.sequence = check_alphabet(sequence, material=material)
        if material is 'dna':
            self.sequence = translate_seq(self.sequence)
        self.frequencies = codon_freq_sc_nested
        self.threshold = threshold

    def __repr__(self):
        return 'RandomCodons generator for %s' % self.sequence

    def generate(self):
        new_table = {}
        for key, value in self.frequencies.iteritems():
            average = 1.0 / len(value)
            vals = value.iteritems()
            new_value = [i for i, x in vals if x / average >= self.threshold]
            if len(new_value) is 0:
                raise ValueError('The threshold has been set so high that it \
                                  excludes all codons of a given amino acid.')
            new_table[key] = new_value
        coding_sequence = [random.choice(new_table[x]) for x in self.sequence]
        coding_sequence = ''.join(coding_sequence)
        return coding_sequence
