'''
Tests for utils submodule of the analysis module.

'''

from nose.tools import assert_equal, assert_raises
from pymbt.analysis import utils
from pymbt import sequence


def test_utils():
    test_DNA = sequence.DNA('ATAGCGATACGAT')
    test_RNA = sequence.RNA('AUGCGAUAGCGAU')
    test_peptide = sequence.Peptide('msvkkkpvqg')
    test_str = 'msvkkkpvgq'

    assert_equal(utils.sequence_type(test_DNA), 'dna')
    assert_equal(utils.sequence_type(test_RNA), 'rna')
    assert_equal(utils.sequence_type(test_peptide), 'peptide')
    assert_raises(Exception, utils.sequence_type, test_str)
