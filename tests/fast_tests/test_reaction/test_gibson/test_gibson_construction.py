import os
from nose.tools import assert_equal
from pymbt import seqio, reaction


def test_construction():
    plasmid_path = os.path.join(os.path.dirname(__file__), 'gibson_test.fasta')
    f1_path = os.path.join(os.path.dirname(__file__), 'fragment_1.fasta')
    f2_path = os.path.join(os.path.dirname(__file__), 'fragment_2.fasta')
    f3_path = os.path.join(os.path.dirname(__file__), 'fragment_3.fasta')
    plasmid = seqio.read_dna(plasmid_path)
    f1 = seqio.read_dna(f1_path)
    f2 = seqio.read_dna(f2_path)
    f3 = seqio.read_dna(f3_path)

    gibsoned = reaction.Gibson([f1, f2, f3]).run()

    expected_length = len(plasmid)
    gibsoned_length = len(gibsoned)
    assert_equal(gibsoned_length, expected_length)
