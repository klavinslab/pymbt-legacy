'''Generate a random DNA sequence.'''

import random
from pymbt import data
from pymbt import sequence


def random_dna(length):
    '''Generate a random DNA sequence.

    :param length: Output sequence length.
    :type length: int

    '''
    return sequence.DNA(''.join(random.choice('ATGC') for i in range(length)))


def random_codons(peptide, frequency_cutoff=0.0, weighted=False,
                  table=data.common.CODON_FREQ_BY_AA['sc']):
    '''Generate randomized codons given a peptide sequence.

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
    # Process codon table using frequency_cutoff
    new_table = {}
    # IDEA: cutoff should be relative to most-frequent codon, not average?
    for aa, codons in table.iteritems():
        average = sum(codons.values()) / len(codons)
        new_table[aa] = {}
        for codon, frequency in codons.iteritems():
            if frequency > frequency_cutoff * average:
                new_table[aa][codon] = frequency

    rna = ''
    for amino_acid in str(peptide):
        codons = new_table[amino_acid.upper()]
        if not codons:
            msg = 'Frequency cutoff prevents use of {}'.format(amino_acid)
            raise ValueError(msg)
        if weighted:
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
        rna += selection
    return sequence.RNA(rna)
