'''
Tests for RandomBases class of analysis module.

'''

from nose.tools import assert_equal
from pymbt import design


def test_randomdna():
    '''
    This test is pretty basic right now - not sure how much checking
    can be done for a random DNA base generator.

    '''

    output = design.RandomDNA(200).run()
    assert_equal(len(output), 200)
