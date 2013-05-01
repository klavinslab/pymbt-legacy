'''Tests primer design module.'''

from nose.tools import assert_equals
from pymbt.primer_design import design_primer
from pymbt.primer_design import design_primer_gene


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
    primer1 = design_primer(seq, tm=72, min_len=10, tm_undershoot=1,
                            tm_overshoot=3, end_gc=False,
                            tm_method='finnzymes', overhang='')
    assert_equals(primer1[0], 'atggtgagcaagggcgaggag')


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
    primers = design_primer_gene(seq, tm=72, min_len=10, tm_undershoot=1,
                                 tm_overshoot=3, end_gc=False,
                                 tm_method='finnzymes', overhangs='')
    assert_equals(primers, [('atggtgagcaagggcgaggag', 72.0133923526509),
                  ('ttacttgtacagctcgtccatgccg', 71.69372970141251)])
