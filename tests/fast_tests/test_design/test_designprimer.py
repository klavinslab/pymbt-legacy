'''Tests primer design module.'''

from nose.tools import assert_equals
from pymbt import design
from pymbt import sequence


def test_design_primer():
    '''Test design_primer function.'''

    seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggc' + \
          'gacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaag' + \
          'ctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgacc' + \
          'accttcggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttc' + \
          'ttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggc' + \
          'aactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctg' + \
          'aagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaac' + \
          'agccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatc' + \
          'cgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatc' + \
          'ggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaa' + \
          'gaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcact' + \
          'ctcggcatggacgagctgtacaagtaa'
    dna_seq = sequence.DNA(seq)
    primer, tm = design.DesignPrimer(dna_seq, tm=72, min_len=10,
                                     tm_undershoot=1, tm_overshoot=3,
                                     end_gc=False, tm_method='finnzymes',
                                     overhang='').run()
    assert_equals(str(primer), 'atggtgagcaagggcgaggag')


def test_design_primer_gene():
    '''Test design_primer_gene function.'''

    seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggc' + \
          'gacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaag' + \
          'ctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgacc' + \
          'accttcggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttc' + \
          'ttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggc' + \
          'aactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctg' + \
          'aagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaac' + \
          'agccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatc' + \
          'cgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatc' + \
          'ggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaa' + \
          'gaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcact' + \
          'ctcggcatggacgagctgtacaagtaa'
    dna_seq = sequence.DNA(seq)
    primers_list = design.DesignPrimerGene(dna_seq, tm=72, min_len=10,
                                           tm_undershoot=1, tm_overshoot=3,
                                           end_gc=False, tm_method='finnzymes',
                                           overhangs='').run()
    primers = [str(x[0]) for x in primers_list]
    assert_equals(primers, ['atggtgagcaagggcgaggag',
                            'ttacttgtacagctcgtccatgccg'])
