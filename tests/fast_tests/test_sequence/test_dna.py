'''
Tests for the DNA sequence class.

'''

from pymbt import sequence
from nose.tools import assert_equal


def test_complements():
    test_dna = sequence.DNA('atgc')

    assert_equal(test_dna.reverse_complement().top, 'gcat')
