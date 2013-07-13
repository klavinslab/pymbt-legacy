'''
Tests for Needleman-Wunsch function 'needle'

'''

from nose.tools import assert_equal
from pymbt.analysis._sequencing.needle import needle


def test_needle():
    ref_seq = 'ATGCGATACGATA'
    res_seqs = ['GATCGATATCGATAT', 'GATACGATATGACGATA', 'TAATTACGGAT']
    expected_alignments = [('-ATGCGATA-CGATA-', 'GAT-CGATATCGATAT'),
                           ('-ATGCGAT---ACGATA', 'GATACGATATGACGATA'),
                           ('--ATGCGATAC-GATA', 'TAAT----TACGGAT-')]
    expected_scores = [40.0, 45.0, 18.5]

    for i, seq in enumerate(res_seqs):
        ref, res, score = needle(ref_seq, seq, gapopen=10, gapextend=0.5)
        assert_equal((ref, res), expected_alignments[i])
        assert_equal(score, expected_scores[i])
