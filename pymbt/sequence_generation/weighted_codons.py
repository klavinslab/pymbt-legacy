import random
from pymbt.common_data import codon_table
from pymbt.common_data import codon_freq
from pymbt.sequence_manipulation import check_alphabet
from pymbt.sequence_manipulation import translate_seq


class WeightedCodons:
    '''Provides a generator class for random DNA or RNA sequence.
         sequence: sequence for which to generate randomized codons.
         frequency_table: specifies which codon frequency dictionary to use.
         material: input material.
             'dna' for DNA or
             'pep' for amino acid sequence (default).'''
    def __init__(self, sequence, frequency_table='sc', material='pep'):
        self.sequence = check_alphabet(sequence, material=material)
        if material == 'dna':
            self.sequence = translate_seq(self.sequence)
        self.material = material
        self.codons = codon_table
        self.codon_freq = codon_freq[frequency_table]

    def __repr__(self):
        return 'RandomCodons generator for %s' % self.sequence

    def weighted(self, pep):
        # Takes an amino acid and selects one at random, weighted by frequency
        cur_codons = self.codons[pep]
        cur_freqs = [self.codon_freq[x] for x in cur_codons]
        cumsum = []
        c = 0
        for i, x in enumerate(cur_freqs):
            c += x
            cumsum.append(c)
        # Using max val instead of 1 - might sum to slightly less than 1
        random_num = random.uniform(0, max(cumsum))
        for i, v in enumerate(cumsum):
            if v > random_num:
                return cur_codons[i]

    def generate(self):
        coding_sequence = [self.weighted(x) for x in self.sequence]
        coding_sequence = ''.join(coding_sequence)

        return coding_sequence
