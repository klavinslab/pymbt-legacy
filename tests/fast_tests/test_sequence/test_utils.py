'''
Tests utils submodule of sequence module.

'''

from nose.tools import assert_equals
from pymbt import sequence


def test_reverse_complement():
    '''Tests reverse_complement function.'''

    seq = 'ATGCNatgcn'
    rev_seq = 'ngcatNGCAT'
    assert_equals(sequence.utils.reverse_complement(seq, 'dna'), rev_seq)
