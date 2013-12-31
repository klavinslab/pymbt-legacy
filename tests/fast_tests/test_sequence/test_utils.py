'''
Tests utils submodule of sequence module.

'''

from nose.tools import assert_equal, assert_raises
from pymbt._sequence import utils


def test_reverse_complement():
    '''
    Tests reverse_complement function.

    '''

    seq = 'ATGCNatgcn'
    rev_seq = 'ngcatNGCAT'
    assert_equal(utils.reverse_complement(seq, 'dna'), rev_seq)


def test_check_alphabet():
    '''
    Tests alphabet checking function.

    '''

    seq = 'ATGCNatgcn'
    bad_seq = 'duck'

    try:
        utils.check_alphabet(seq, 'dna')
    except:
        # HACK
        assert False
    assert_raises(ValueError, utils.check_alphabet, seq, 'derp')
    assert_raises(ValueError, utils.check_alphabet, bad_seq, 'dna')
    assert_raises(ValueError, utils.check_alphabet, bad_seq, 'rna')
    assert_raises(ValueError, utils.check_alphabet, bad_seq, 'peptide')
