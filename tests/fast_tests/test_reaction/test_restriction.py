'''Test restriction reaction module.'''

from nose.tools import assert_equal
from pymbt import sequence
from pymbt import reaction


class TestDigest(object):
    '''Test digest function.'''
    def __init__(self):
        # Contains NcoI site
        self.dna = sequence.DNA('TGACCATGGAAA')

    def test_not_found(self):
        '''If site not found, should return input sequence in list.'''
        ecorv = sequence.RestrictionSite(sequence.DNA('GATATC'), (3, 3),
                                         name='EcoRV')
        assert_equal(self.dna, reaction.digest(self.dna, ecorv)[0])

    def test_ncoi_cut(self):
        '''Test standard TypeII cutter.'''
        ncoi = sequence.RestrictionSite(sequence.DNA('CCATGG'), (1, 5),
                                        name='NcoI')
        assert_equal(reaction.digest(self.dna, ncoi),
                     [sequence.DNA('TGAC----', bottom='CATGGTCA'),
                      sequence.DNA('CATGGAAA', bottom='TTTC----')])
        assert_equal(reaction.digest(self.dna.circularize(), ncoi),
                     [sequence.DNA('CATGGAAATGAC----',
                                   bottom='CATGGTCATTTC----')])

    def test_ecorv_cut(self):
        '''Test blunt-end cutter.'''
        ecorv = sequence.RestrictionSite(sequence.DNA('GATATC'), (3, 3),
                                         name='EcoRV')
        assert_equal(reaction.digest(sequence.DNA('GATATC'), ecorv),
                     [sequence.DNA('GAT'), sequence.DNA('ATC')])

    def test_psti_cut(self):
        '''Test 3\' cutter.'''
        psti = sequence.RestrictionSite(sequence.DNA('CTGCAG'), (5, 1),
                                        name='PstI')
        assert_equal(reaction.digest(sequence.DNA('ACTGCAGA'), psti),
                     [sequence.DNA('ACTGCA', bottom='----GT'),
                      sequence.DNA('----GA', bottom='TCTGCA')])
