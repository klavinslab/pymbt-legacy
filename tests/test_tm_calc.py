'''Tests melting temp module (tm_calc).'''

from nose.tools import assert_equals
from pymbt.analysis import Tm
from pymbt.sequence.dna import DNA


def test_finnzymes():
    '''Tests finnzymes method output.'''

    melt = Tm(DNA('ATGCGATAGCGATAGC'), method='finnzymes').run()
    assert_equals(melt, 55.237003002075255)
