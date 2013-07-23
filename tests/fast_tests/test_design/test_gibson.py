'''Test gibson design module.'''
from nose.tools import assert_equal
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
    # Design primers - with homology all on left side, right side, or mixed
    left = design.gibson_primers(tdh3_3prime, yfp_nterm, 'left')
    left_primer1 = sequence.Primer(sequence.DNA('gaagtttaagaaatagatttacagaat' +
                                                'tacaatcaatacctaccgtcttt'),
                                   tm=71.99388428631767)
    left_primer2 = sequence.Primer(sequence.DNA('atggtgagcaagggcgaggag'),
                                   tm=72.01339235265078,
                                   overhang=sequence.DNA('GGTATTGATTGTAATTCT' +
                                                         'GTAAATCTATTTCTTAAA' +
                                                         'CTTC'))
    assert_equal(left, (left_primer1, left_primer2))

    right = design.gibson_primers(tdh3_3prime, yfp_nterm, 'right')
    right_primer1 = sequence.Primer(sequence.DNA('gaagtttaagaaatagatttacagaa' +
                                                 'ttacaatcaatacctaccgtcttt'),
                                    tm=71.99388428631767,
                                    overhang=sequence.DNA('TCGCCCTTGCTCACCAT'))
    right_primer2 = sequence.Primer(sequence.DNA('atggtgagcaagggcgaggag'),
                                    tm=72.01339235265078)
    assert_equal(right, (right_primer1, right_primer2))

    mixed = design.gibson_primers(tdh3_3prime, yfp_nterm, 'mixed')
    mixed_primer1 = sequence.Primer(sequence.DNA('gaagtttaagaaatagatttacagaa' +
                                                 'ttacaatcaatacctaccgtcttt'),
                                    tm=71.99388428631767,
                                    overhang=sequence.DNA('ATGGTGAGCAAGG'))
    mixed_primer2 = sequence.Primer(sequence.DNA('atggtgagcaagggcgaggag'),
                                    tm=72.01339235265078,
                                    overhang=sequence.DNA('TTCTTAAACTTC'))
    assert_equal(mixed, (mixed_primer1, mixed_primer2))
