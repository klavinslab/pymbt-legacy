'''
Generate a random DNA sequence.

'''

import random
from pymbt import data
from pymbt import sequence


class RandomDNA(object):
    '''
    Generate a random DNA sequence.

    '''

    def __init__(self, length):
        '''
        :param length: Output sequence length.
        :type length: int

        '''

        self.length = length

    def __repr__(self):
        '''
        Printed representation.

        '''

        return '{} of length {}.'.format(self.__class__.__name__, self.length)

    def random_base(self):
        '''
        Generate a random base.

        '''

        return sequence.DNA(random.choice('ATGC'))

    def run(self):
        '''
        Generate the sequence.

        '''

        result = sum([self.random_base() for i in range(self.length)])
        return result


class RandomCodons(object):
    '''
    Generator class for random DNA or RNA sequence.

    '''

    def __init__(self, peptide, frequency_cutoff=0.0, weighted=False,
                 table=data.common.CODON_FREQ_BY_AA['sc']):
        '''
        :param peptide: Peptide sequence for which to generate randomized
                        codons.
        :type peptide: pymbt.sequence.Peptide
        :param frequency_cutoff: Relative codon usage cutoff - codons that
                                 are rarer will not be used. Frequency is
                                 relative to average over all codons for a
                                 given amino acid.
        :param frequency_cutoff: Codon frequency table to use.
        :param weighted: Use codon table
        :type weighted: bool
        :param table: Codon frequency table to use. Table should be organized
                      by amino acid, then be a dict of codon: frequency.
                      Only relevant if weighted=True or frequency_cutoff > 0.
                      Tables available:

                      data.CODON_FREQ_SC_NESTED
        :type table: dict

        '''

        self.template = peptide
        self.frequency_cutoff = frequency_cutoff
        self.weighted = weighted
        # Process table given frequency_cutoff
        self.table = {}
        # IDEA: cutoff should be relative to most-frequent codon, not average?
        for aa, codons in table.iteritems():
            average = sum(codons.values()) / len(codons)
            self.table[aa] = {}
            for codon, frequency in codons.iteritems():
                if frequency > frequency_cutoff * average:
                    self.table[aa][codon] = frequency

    def __repr__(self):
        '''
        Printed representation.

        '''

        return '{0} for {1}'.format(self.__class__.__name__, self.template)

    def choose_codon(self, amino_acid):
        codons = self.table[amino_acid.upper()]
        if not codons:
            msg = 'Frequency cutoff prevents use of {}'.format(amino_acid)
            raise ValueError(msg)
        if self.weighted:
            cumsum = []
            running_sum = 0
            for codon, frequency in codons.iteritems():
                running_sum += frequency
                cumsum.append(running_sum)
            # Using max val instead of 1 - might sum to slightly less than 1
            random_num = random.uniform(0, max(cumsum))
            for codon, value in zip(codons, cumsum):
                if value > random_num:
                    selection = codon
                    break
        else:
            selection = random.choice(codons.keys())

        return sequence.RNA(selection)

    def run(self):
        '''
        Generate the sequence.

        '''

        result = sum([self.choose_codon(str(a)) for a in self.template])
        return result
