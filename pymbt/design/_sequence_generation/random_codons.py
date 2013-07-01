'''
Generate random codons.

'''

import random
from pymbt.data.common import CODON_FREQ_SC_NESTED
from pymbt import reaction
from pymbt import sequence


class RandomCodons(object):
    '''
    Generator class for random DNA or RNA sequence.

    '''

    def __init__(self, dna_object, threshold=0.0):
        '''
        :param sequence: Sequence for which to generate randomized codons.
        :type sequence: str
        :param freq_table: Codon frequency table to use.
        :type freq_table: dict
        :param threshold: relative threshold cutoff for codon frequencies.
        :type threshold: float between 0 and 1

        '''

        self.template = dna_object
        self.rna = reaction.Transcription(self.template).run()
        self.peptide = reaction.Translation(self.rna).run()
        self.frequencies = CODON_FREQ_SC_NESTED
        self.threshold = threshold

    def __repr__(self):
        '''
        Printed representation of RandomCodons object.

        '''
        return 'RandomCodons generator for {}'.format(self.peptide)

    def run(self):
        '''
        Generate the sequence.

        '''

        new_table = {}
        for key, value in self.frequencies.iteritems():
            average = 1.0 / len(value)
            vals = value.iteritems()
            new_value = [i for i, x in vals if x / average >= self.threshold]
            if not new_value:
                raise ValueError('The threshold has been set so high that it \
                                  excludes all codons of a given amino acid.')
            new_table[key] = new_value
        coding_sequence = [random.choice(new_table[str(x).upper()]) for x in
                           self.peptide]
        coding_rna = sequence.RNA(''.join(coding_sequence))
        coding_dna = reaction.ReverseTranscription(coding_rna).run()
        return coding_dna
