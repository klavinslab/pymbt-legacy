'''
Tests for RandomCodons class of analysis module.

'''

from nose.tools import assert_equal
from nose.tools import assert_not_equal
from pymbt import reaction
from pymbt import design


def test_randomcodons():
    '''
    This test is pretty basic right now - not sure how much checking
    can be done for a random DNA base generator.

    '''

    reference_seq = 'AUGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAG'
    reference_peptide = reaction.Translation(reference_seq).run()
    output = design.RandomCodons(reference_peptide).run()
    output_peptide = reaction.Translation(reference_seq).run()

    assert_equal(len(output), len(reference_seq))
    assert_equal(reference_peptide, output_peptide)
    assert_not_equal(reference_seq, output)
