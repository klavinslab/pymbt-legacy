# TODO:
# Check that the input sequence is amino acids
# Haven't yet implemented 'dna' functionality

import random
import pymbt.sequence_generation.codon_data


class WeightedCodons:
    '''Provides a generator class for random DNA or RNA sequence.
         sequence: sequence for which to generate randomized codons.
         material: 'dna' for DNA or 'aa' for amino acid sequence.
    '''
    def __init__(self, sequence, organism='sc', material='aa'):
        self.sequence = sequence
        self.material = material
        self.codons = codon_data.codontable
        self.codon_freq = codon_data.codon_freq[organism]
        # cumulative sums of amino acids - useful for weights

    def __repr__(self):
        return 'RandomCodons generator for %s' % self.sequence

    def weighted(self, aa):
        # Takes an amino acid and selects one at random, weighted by frequency
        cur_codons = self.codons[aa]
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
        if self.material is not 'dna' and self.material is not 'aa':
            raise ValueError("material must be 'dna' or 'aa")
        coding_sequence = []
        coding_sequence = [self.weighted(x) for x in self.sequence]
        coding_sequence = ''.join(coding_sequence)
        return coding_sequence
