'''Test gibson design module.'''
from nose.tools import assert_equal, assert_raises
from pymbt import design, sequence


def test_gibson_primers():
    '''Test gibson_primers function.'''
    # Fuse tdh3 promoter sequence to yfp (trimmed for readability)
    tdh3_3prime = sequence.DNA('aaccagttccctgaaattattcccctacttgactaataagtat' +
                               'ataaagacggtaggtattgattgtaattctgtaaatctatttc' +
                               'ttaaacttc')
    yfp_nterm = sequence.DNA('atggtgagcaagggcgaggagctgttcaccggggtggtgcccatc' +
                             'ctggtcgagctggacggcgacgtaaacggccacaagttcagcgtg' +
                             'tccggcgagggcgagggcgatgccacctacggcaagctgaccctg' +
                             'aag')
    # Expected annealing sequences and their Tms
    fwd_anneal = sequence.DNA('atggtgagcaagggcg')
    fwd_tm = 64.64172107821065
    rev_anneal = sequence.DNA('gaagtttaagaaatagatttacagaattacaatcaatac')
    rev_tm = 64.24536287254085
    # Expected overlaps
    all_right = sequence.DNA('TCGCCCTTGCTCACCAT')
    all_left = sequence.DNA('GGTATTGATTGTAATTCTGTAAATCTATTTCTTAAACTTC')
    mixed_fwd = sequence.DNA('TTCTTAAACTTC')
    mixed_rev = sequence.DNA('CCTTGCTCACCAT')
    # Design primers - with homology all on left side, right side, or mixed
    # All on the 'right' - i.e. fwd primer
    right = design.gibson_primers(tdh3_3prime, yfp_nterm, 'right')
    right_rev = sequence.Primer(rev_anneal, tm=rev_tm, overhang=all_right)
    right_fwd = sequence.Primer(fwd_anneal, tm=fwd_tm)
    assert_equal(right, (right_rev, right_fwd))
    # All on the 'left' - i.e. rev primer
    left = design.gibson_primers(tdh3_3prime, yfp_nterm, 'left')
    left_rev = sequence.Primer(rev_anneal, tm=rev_tm)
    left_fwd = sequence.Primer(fwd_anneal, tm=fwd_tm, overhang=all_left)
    assert_equal(left, (left_rev, left_fwd))
    # On both primers
    mixed = design.gibson_primers(tdh3_3prime, yfp_nterm, 'mixed')
    mixed_primer1 = sequence.Primer(rev_anneal, tm=rev_tm, overhang=mixed_rev)
    mixed_primer2 = sequence.Primer(fwd_anneal, tm=fwd_tm, overhang=mixed_fwd)
    assert_equal(mixed, (mixed_primer1, mixed_primer2))

    assert_raises(ValueError, design.gibson_primers, tdh3_3prime, yfp_nterm,
                  'duck')
