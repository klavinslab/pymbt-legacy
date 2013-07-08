'''
Tests for the Tm analysis class.

'''

from nose.tools import assert_equal
from pymbt import analysis
from pymbt import sequence


def test_finnzymes():
    '''
    Tests finnzymes method output.

    '''

    melt = analysis.tm(sequence.DNA('ATGCGATAGCGATAGC'), parameters='cloning')
    assert_equal(melt, 55.2370030020752)
