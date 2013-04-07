import random
from pymbt.common_data import codon_freq_sc_nested
from pymbt.sequence_manipulation import check_alphabet
from pymbt.sequence_manipulation import translate_seq


class RandomCodons:
    '''Generator class for random DNA or RNA sequence.'''
    def __init__(self, sequence, material='pep', freq_table='sc',
                 threshold=0.0):
        '''
        :param sequence: Sequence for which to generate randomized codons.
        :type sequence: str.
        :param material: 'dna' for DNA or 'pep' for amino acid sequence.
        :type material: str.

        '''
        self.sequence = check_alphabet(sequence, material=material)
        if material == 'dna':
            self.sequence = translate_seq(self.sequence)
        self.frequencies = codon_freq_sc_nested
        self.threshold = threshold

    def __repr__(self):
        return 'RandomCodons generator for %s' % self.sequence

    def generate(self):
        '''Generate the sequence.'''
        new_table = {}
        for key, value in self.frequencies.iteritems():
            average = 1.0 / len(value)
            vals = value.iteritems()
            new_value = [i for i, x in vals if x / average >= self.threshold]
            if not new_value:
                raise ValueError('The threshold has been set so high that it \
                                  excludes all codons of a given amino acid.')
            new_table[key] = new_value
        coding_sequence = [random.choice(new_table[x]) for x in self.sequence]
        coding_sequence = ''.join(coding_sequence)
        return coding_sequence
