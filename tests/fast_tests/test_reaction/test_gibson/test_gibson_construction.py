import os
from nose.tools import assert_equal
from pymbt import seqio, reaction


def test_construction():
    plasmid_path = os.path.join(os.path.dirname(__file__), 'gibson_test.fasta')
    f1_path = os.path.join(os.path.dirname(__file__), 'fragment_1.fasta')
    f2_path = os.path.join(os.path.dirname(__file__), 'fragment_2.fasta')
    f3_path = os.path.join(os.path.dirname(__file__), 'fragment_3.fasta')
    f3_linear_path = os.path.join(os.path.dirname(__file__),
                                  'fragment_3_linear.fasta')
    plasmid = seqio.read_dna(plasmid_path)
    f1 = seqio.read_dna(f1_path)
    f2 = seqio.read_dna(f2_path)
    f3 = seqio.read_dna(f3_path)
    f3_linear = seqio.read_dna(f3_linear_path)

    gibsoned_circular = reaction.Gibson([f1, f2, f3]).run()
    gibsoned_linear = reaction.Gibson([f1, f2, f3_linear]).run(linear=True)

    expected_length = len(plasmid)
    gibsoned_circular_length = len(gibsoned_circular)
    gibsoned_linear_length = len(gibsoned_linear)
    assert_equal(gibsoned_circular_length, expected_length)
    assert_equal(gibsoned_linear_length, expected_length)
    assert_equal(gibsoned_circular.topology, 'circular')
    assert_equal(gibsoned_linear.topology, 'linear')
    assert_equal(plasmid.top, gibsoned_circular.top)
    assert_equal(plasmid.top, gibsoned_linear.top)
