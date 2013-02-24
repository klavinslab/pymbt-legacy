# TODO:
# Check that the input sequence is amino acids
# Haven't yet implemented 'dna' functionality
# Threshold number is a bit strange - a value of '1.0' excludes
# only codons that fall below average
# Put thresholding function in a separate file

import random
import pymbt.sequence_generation.codon_data


class RandomCodons:
    '''Provides a generator class for random DNA or RNA sequence.
         sequence: sequence for which to generate randomized codons.
         material: 'dna' for DNA or 'aa' for amino acid sequence.
    '''
    def __init__(self, sequence, organism, material='aa', threshold=0.0):
        self.sequence = sequence
        self.material = material
        self.organism = organism
        self.threshold = threshold
        self.codons = codon_data.codontable
        self.codon_freqs = codon_data.codon_freq[self.organism]

    def __repr__(self):
        return 'RandomCodons generator for %s' % self.sequence

    def generate(self):
        if self.material is not 'dna' and self.material is not 'aa':
            raise ValueError("material must be 'dna' or 'aa")
        mod_codons = self.codons.copy()
        for key, val in mod_codons.iteritems():
            freq_list = self.codon_freqs
            for x in val:
                c_freq = self.codon_freqs[x]
                threshold = self.threshold / len(val)
                mod_codons[key] = [x for x in val if c_freq >= threshold]
            if len(mod_codons[key]) is 0:
                raise ValueError('The threshold has been set so high that it
                                  excludes all codons of a given amino acid.')
        coding_sequence = [random.choice(mod_codons[x]) for x in self.sequence]
        coding_sequence = ''.join(coding_sequence)
        return coding_sequence
