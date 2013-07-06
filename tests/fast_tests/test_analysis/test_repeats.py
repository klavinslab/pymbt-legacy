'''
Tests Repeats analysis class.

'''

from nose.tools import assert_equal
from pymbt import analysis
from pymbt import sequence


def test_find_repeats():
    reference_seq = sequence.DNA('atgatgccccgatagtagtagtag')
    reference_result = [('atg', 2), ('gat', 2), ('tag', 4), ('gta', 3),
                        ('agt', 3), ('ccc', 2)]

    output_result = analysis.Repeats(reference_seq, 3).run()
    assert_equal(output_result, reference_result)
