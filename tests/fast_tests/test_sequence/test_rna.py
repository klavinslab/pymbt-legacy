'''
Tests for the RNA sequence class.

'''

from pymbt import sequence
from nose.tools import assert_equal


def test_complements():
    test_rna = sequence.RNA('augc')

    assert_equal(test_rna.reverse_complement().top, 'gcau')
