import os
from nose.tools import assert_equal, assert_raises
from pymbt import reaction, seqio


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

    gibsoned_circular = reaction.Gibson([f1, f2, f3]).run_circular()
    gibsoned_linear = reaction.Gibson([f1, f2, f3_linear]).run_linear()

    expected_length = len(plasmid)
    gibsoned_circular_length = len(gibsoned_circular)
    gibsoned_linear_length = len(gibsoned_linear)
    assert_equal(gibsoned_circular_length, expected_length)
    assert_equal(gibsoned_linear_length, expected_length)
    assert_equal(gibsoned_circular.topology, 'circular')
    assert_equal(gibsoned_linear.topology, 'linear')
    assert_equal(str(plasmid), str(gibsoned_circular))
    assert_equal(str(plasmid), str(gibsoned_linear))

    # Should fail with circular input
    assert_raises(ValueError, reaction.Gibson, [f1.circularize()])
    # Should fail if compatible end can't be found
    assert_raises(Exception, reaction.Gibson([f1, f3[30:-30]]).run_linear)
    normal = [f1, f2, f3]
    rotated = [f1, f2, f3.reverse_complement()]
    # Gibson should work regardless of fragment orientation
    assert_equal(reaction.Gibson(normal).run_circular(),
                 reaction.Gibson(rotated).run_circular())
    # A redundant fragment shouldn't affect the outcome
    assert_equal(reaction.Gibson([f1, f2, f3]).run_circular(),
                 reaction.Gibson([f1, f2, f2, f3]).run_circular())
    # A fragment that can't circularize should raise a ValueError
    assert_raises(ValueError, reaction.Gibson([f1, f2, f3[:-80]]).run_circular)
    # But should still work fine as a linear fragment
    assert_equal(reaction.Gibson([f1, f2, f3]).run_linear()[:-80],
                 reaction.Gibson([f1, f2, f3[:-80]]).run_linear())
    # If there's more than one way to make the Gibson happen, should error
    assert_raises(reaction._gibson.AmbiguousGibsonError,
                  reaction.Gibson([f1, f2, f2[:60] + f3, f3]).run_circular)
