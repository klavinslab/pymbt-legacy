from nose.tools import assert_equals
from pymbt import tm_calc


def test_finnzymes():
    melt = tm_calc.calc_tm('ATGCGATAGCGATAGC', method='finnzymes')
    assert_equals(melt, 55.237003002075255)
