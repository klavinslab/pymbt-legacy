'''
Generate random, but usage frequency-weighted codons (i.e. codon optimization).

'''

import random
from pymbt.data.common_data import CODON_TABLE, CODON_FREQ
from pymbt.sequence_utils import translate_seq
from pymbt import sequence


class WeightedCodons(object):
    '''Provides a generator class for random, weighted DNA or RNA sequences.'''

    def __init__(self, dna_object, frequency_table='sc', material='dna'):
        '''
        :param sequence: Input sequence.
        :type sequence: str
        :param frequency_table: The codon frequency dictionary to use.
        :type frequency_table: dict
        :param material: 'dna' for DNA.
        :type material: str

        '''

        self.template = str(dna_object)
        self.pep = translate_seq(self.template)
        self.material = material
        self.codons = CODON_TABLE
        self.codon_freq = CODON_FREQ[frequency_table]

    def __repr__(self):
        return 'RandomCodons generator for {}'.format(self.pep)

    def weighted(self, pep):
        '''
        Take an amino acid, select a codon at random, weighted by frequency.

        :param pep: Peptide sequence.
        :type pep: str

        '''

        codons = self.codons[pep]
        frequencies = [self.codon_freq[x] for x in codons]
        cumsum = []
        running_sum = 0
        for i, frequency in enumerate(frequencies):
            running_sum += frequency
            cumsum.append(running_sum)
        # Using max val instead of 1 - might sum to slightly less than 1
        random_num = random.uniform(0, max(cumsum))
        for i, value in enumerate(cumsum):
            if value > random_num:
                return codons[i]

    def generate(self):
        '''Generate the sequence.'''

        coding_sequence = [self.weighted(x) for x in self.pep]
        coding_sequence = sequence.DNA(''.join(coding_sequence))

        return coding_sequence
