'''
Generate random, but usage frequency-weighted codons (i.e. codon optimization).

'''

import random
from pymbt.data.common import CODON_TABLE, CODON_FREQ
from pymbt import reaction
from pymbt import sequence


class WeightedCodons(object):
    '''
    Provides a generator class for random, weighted DNA or RNA sequences.

    '''

    def __init__(self, dna_object, frequency_table='sc'):
        '''
        :param sequence: Input sequence.
        :type sequence: str
        :param frequency_table: The codon frequency dictionary to use.
        :type frequency_table: dict

        '''

        self.template = dna_object
        self.rna = reaction.Transcription(self.template).run()
        self.peptide = reaction.Translation(self.rna).run()
        self.codons = CODON_TABLE
        self.codon_freq = CODON_FREQ[frequency_table]

    def __repr__(self):
        '''
        Printed representation of WeightedCodons object.

        '''

        return 'RandomCodons generator for {}'.format(self.peptide)

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
        for frequency in frequencies:
            running_sum += frequency
            cumsum.append(running_sum)
        # Using max val instead of 1 - might sum to slightly less than 1
        random_num = random.uniform(0, max(cumsum))
        for codon, value in zip(codons, cumsum):
            if value > random_num:
                return codon

    def run(self):
        '''
        Generate the sequence.

        '''

        coding_sequence = [self.weighted(str(x).upper()) for x in self.peptide]
        coding_rna = sequence.RNA(''.join(coding_sequence))
        coding_dna = reaction.ReverseTranscription(coding_rna).run()

        return coding_dna
