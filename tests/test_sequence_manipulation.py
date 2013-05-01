'''Tests sequence_manipulation module.'''

from nose.tools import assert_equals
from pymbt import sequence_manipulation


def test_reverse_complement():
    '''Tests reverse_complement function.'''

    seq = 'ATGCNatgcn'
    rev_seq = 'ngcatNGCAT'
    assert_equals(sequence_manipulation.reverse_complement(seq), rev_seq)


def test_translate_seq():
    '''Tests DNA translation function.'''

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
    prot = 'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLV' + \
           'TTFGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRI' + \
           'ELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQN' + \
           'TPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK'
    trans = sequence_manipulation.translate_seq(seq)
    assert_equals(trans, prot)
