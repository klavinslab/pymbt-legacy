'''
Tests for Needleman-Wunsch function 'needle'

'''

from nose.tools import assert_equal
from pymbt.analysis._sequencing.needle import needle


def test_needle():
    ref_seq = 'ATGCGATACGATA'
    res_seqs = ['GATCGATATCGATAT', 'GATACGATATGACGATA', 'TAATTACGGAT']

    aligned = needle(ref_seq, res_seqs, gapopen=10, gapextend=0.5)
    alignments = aligned['alignments']
    scores = aligned['scores']

    expected_alignments = [('-ATGCGATA-CGATA-', 'GAT-CGATATCGATAT'),
                           ('-ATGCGAT---ACGATA', 'GATACGATATGACGATA'),
                           ('--ATGCGATAC-GATA', 'TAAT----TACGGAT-')]
    expected_scores = [40.0, 45.0, 18.5]

    assert_equal(alignments, expected_alignments)
    assert_equal(scores, expected_scores)
