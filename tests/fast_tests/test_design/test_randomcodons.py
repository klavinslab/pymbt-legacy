'''
Tests for RandomCodons class of analysis module.

'''

from nose.tools import assert_equal
from nose.tools import assert_not_equal
from pymbt import reaction
from pymbt import sequence
from pymbt import design


def test_randomcodons():
    '''
    This test is pretty basic right now - not sure how much checking
    can be done for a random DNA base generator.

    '''

    reference_seq = sequence.RNA('AUGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAG')
    reference_peptide = reaction.translate(reference_seq)
    output = design.random_codons(reference_peptide)
    output_peptide = reaction.translate(reference_seq)

    assert_equal(len(output), len(reference_seq) - 3)
    assert_equal(reference_peptide, output_peptide)
    assert_not_equal(reference_seq, output)
