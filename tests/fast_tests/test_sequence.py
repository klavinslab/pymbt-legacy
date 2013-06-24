from pymbt import sequence
from nose.tools import assert_equal


def test_complements():
    test_dna = sequence.DNA('atgc')
    test_rna = sequence.RNA('augc')

    assert_equal(test_dna.reverse_complement().top, 'gcat')
    assert_equal(test_rna.reverse_complement().top, 'gcau')
