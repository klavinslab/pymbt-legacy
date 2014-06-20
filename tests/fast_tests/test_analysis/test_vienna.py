'''Tests for Vienna RNA module.'''
from nose.tools import assert_equal
from pymbt import analysis, DNA


class TestVienna(object):
    '''Tests for Vienna class.'''
    def __init__(self):
        # M13R primer sequence
        self.dna = DNA('agcggataacaatttcacacaggaaacagctatgaccatg')
        self.vienna = analysis.Vienna([self.dna])

    def test_mfe(self):
        mfe = self.vienna.mfe(temp=37.0)
        assert_equal(mfe, [-0.8])

    def test_unbound_probability(self):
        pairs = self.vienna.unbound_probability(temp=50.0)
        expected = [[0.8423227650608299,
                     0.7364112598806074,
                     0.7739595904918544,
                     0.8080405925870966,
                     0.7687890414629115,
                     0.8313820227858698,
                     0.8841212611807588,
                     0.9171987107347552,
                     0.9940845356599891,
                     0.876512878693969,
                     0.995040828244263,
                     0.9842400552110041,
                     0.7193295801979633,
                     0.6421101022790772,
                     0.5444699866930739,
                     0.5505514660818621,
                     0.9725128307271304,
                     0.9030292003286778,
                     0.9567142431512005,
                     0.9589035939638103,
                     0.9141180723138675,
                     0.7325044907419416,
                     0.5841763402337967,
                     0.6747545663707104,
                     0.6947806160791583,
                     0.7902268103252766,
                     0.8443502137246074,
                     0.9184735503231636,
                     0.702798206583925,
                     0.7065186512835412,
                     0.6989258629215821,
                     0.9195052599358309,
                     0.7655574180108721,
                     0.823148718657445,
                     0.9710836005464278,
                     0.832666101753772,
                     0.8417148542407485,
                     0.9963905980186613,
                     0.9586017081846266,
                     0.9726247543949348]]
        assert_equal(pairs, expected)
