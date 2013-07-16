'''
Tests for Needleman-Wunsch function 'needle'

'''

from nose.tools import assert_equal
from pymbt.analysis._sequencing.needle import needle
from pymbt.sequence import DNA


def test_needle():
    ref_seq = DNA('ATGCGATACGATA')
    res_seqs = [DNA('GATCGATATCGATAT'), DNA('GATACGATATGACGATA'),
                DNA('TAATTACGGAT')]
    expected_alignments = [(DNA('-ATGCGATA-CGATA-'), DNA('GAT-CGATATCGATAT')),
                           (DNA('-ATGCGAT---ACGATA'),
                            DNA('GATACGATATGACGATA')),
                           (DNA('--ATGCGATAC-GATA'), DNA('TAAT----TACGGAT-'))]
    expected_scores = [40.0, 45.0, 18.5]

    for i, seq in enumerate(res_seqs):
        ref, res, score = needle(ref_seq, seq, gapopen=10, gapextend=0.5)
        assert_equal((ref, res), expected_alignments[i])
        assert_equal(score, expected_scores[i])
