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
    primer = design.design_primer(dna_seq, tm=72, min_len=10,
                                  tm_undershoot=1, tm_overshoot=3,
                                  end_gc=False, tm_parameters='cloning',
                                  overhang='')
    assert_equals(str(primer), 'atggtgagcaagggcgaggag')


def test_design_primer_pcr():
    '''Test design_primer_pcr function.'''

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
    primers_list = design.design_primer_pcr(dna_seq, tm=72, min_len=10,
                                            tm_undershoot=1, tm_overshoot=3,
                                            end_gc=False,
                                            tm_parameters='cloning',
                                            overhangs='')
    primers = [str(x.primer()) for x in primers_list]
    assert_equals(primers, ['atggtgagcaagggcgaggag',
                            'ttacttgtacagctcgtccatgccg'])
