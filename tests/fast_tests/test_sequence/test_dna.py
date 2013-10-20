'''Tests for the DNA sequence class.'''
from pymbt import sequence
from pymbt import reaction
from nose.tools import assert_equal, assert_false, assert_true, assert_raises
from nose.tools import assert_not_equal


class TestDNA(object):
    '''Testing class for sequence.DNA'''
    def __init__(self):
        self.test_dna = sequence.DNA('atgc')

    def test_reverse_complement(self):
        assert_equal(str(self.test_dna.reverse_complement()), 'gcat')

    def test_linearize(self):
        circ_dna = self.test_dna.circularize()
        assert_equal(str(circ_dna.linearize()), 'atgc')
        assert_equal(str(circ_dna.linearize(0)), 'atgc')
        assert_equal(str(circ_dna.linearize(1)), 'tgca')
        assert_equal(str(circ_dna.linearize(2)), 'gcat')
        assert_equal(str(circ_dna.linearize(3)), 'catg')
        assert_equal(str(circ_dna.linearize(-1)), 'catg')
        assert_raises(ValueError, circ_dna.linearize().linearize)

    def test_set_stranded(self):
        assert_equal(self.test_dna.set_stranded('ds'), self.test_dna)
        ss_dna = self.test_dna.set_stranded('ss')
        assert_equal(ss_dna.stranded, 'ss')
        ds_dna = self.test_dna.set_stranded('ds')
        assert_equal(ds_dna.stranded, 'ds')
        assert_equal(ds_dna.top(), str(ss_dna))

        ds_to_ss_to_ds = self.test_dna.set_stranded('ss').set_stranded('ds')
        assert_equal(self.test_dna, ds_to_ss_to_ds)

        empty_top = reaction.three_resect(self.test_dna, 400)
        assert_equal(empty_top.set_stranded('ds'), self.test_dna)

        assert_raises(ValueError, self.test_dna.set_stranded, 'duck')

    def test_locate(self):
        assert_equal(self.test_dna.locate('a'), ([0], [2]))
        assert_equal(self.test_dna.locate('at'), ([0], [2]))
        assert_equal(self.test_dna.locate('gc'), ([2], [0]))
        assert_equal(self.test_dna.locate('atgg'), ([], []))
        # Circular DNA tests
        assert_equal(self.test_dna.circularize().locate('a'), ([0], [2]))
        assert_equal(self.test_dna.circularize().locate('at'), ([0], [2]))
        assert_equal(self.test_dna.circularize().locate('gc'), ([2], [0]))
        assert_equal(self.test_dna.circularize().locate('atgg'), ([], []))

    def test_copy(self):
        assert_equal(self.test_dna, self.test_dna.copy())

    def test_palindrome(self):
        palindromic_seq_even = sequence.DNA('ATGCGCAT')
        nonpalindromic_seq_even = sequence.DNA('ATGCGCAA')
        almost_palindrome_odd = sequence.DNA('ATGCCAT')
        assert_true(palindromic_seq_even.is_palindrome())
        assert_false(nonpalindromic_seq_even.is_palindrome())
        assert_false(almost_palindrome_odd.is_palindrome())

    def test_getitem(self):
        assert_equal(str(self.test_dna[0]), 'a')
        assert_equal(str(self.test_dna[1]), 't')
        assert_equal(str(self.test_dna[2]), 'g')
        assert_equal(str(self.test_dna[3]), 'c')
        assert_equal(str(self.test_dna[-1]), 'c')

    def test_delitem(self):
        copy0 = self.test_dna.copy()
        del copy0[0]
        assert_equal(str(copy0), 'tgc')
        copy1 = self.test_dna.copy()
        del copy1[1]
        assert_equal(str(copy1), 'agc')
        copy2 = self.test_dna.copy()
        del copy2[2]
        assert_equal(str(copy2), 'atc')
        copy3 = self.test_dna.copy()
        del copy3[3]
        assert_equal(str(copy3), 'atg')
        copy_1 = self.test_dna.copy()
        del copy_1[-1]
        assert_equal(str(copy_1), 'atg')

    def test_setitem(self):
        copy0 = self.test_dna.copy()
        copy0[0] = 't'
        assert_equal(str(copy0), 'ttgc')
        copy1 = self.test_dna.copy()
        copy1[1] = 'a'
        assert_equal(str(copy1), 'aagc')
        copy2 = self.test_dna.copy()
        copy2[2] = 'a'
        assert_equal(str(copy2), 'atac')
        copy3 = self.test_dna.copy()
        copy3[3] = 'a'
        assert_equal(str(copy3), 'atga')
        copy_1 = self.test_dna.copy()
        copy_1[-1] = 'a'
        assert_equal(str(copy_1), 'atga')

        def set_gap(seq):
            seq[2] = '-'
        assert_raises(ValueError, set_gap, self.test_dna)

    def test_repr(self):
        expected_repr = 'linear dsDNA:\natgc\ntacg'
        assert_equal(repr(self.test_dna), expected_repr)

        expected_circ_repr = 'circular dsDNA:\natgc\ntacg'
        assert_equal(repr(self.test_dna.circularize()), expected_circ_repr)

        repr_1 = 'linear dsDNA:\natgcatgcatgcatgcatgcatgcatgcatgcatgcatgc ... '
        repr_2 = 'atgcatgcatgcatgcatgcatgcatgcatgcatgcatgc\n'
        repr_3 = 'tacgtacgtacgtacgtacgtacgtacgtacgtacgtacg ... '
        repr_4 = 'tacgtacgtacgtacgtacgtacgtacgtacgtacgtacg'
        expected_long_repr = repr_1 + repr_2 + repr_3 + repr_4
        assert_equal(repr(self.test_dna * 50), expected_long_repr)

    def test_str(self):
        assert_equal(str(self.test_dna), 'atgc')

    def test_len(self):
        assert_equal(len(self.test_dna), 4)

    def test_add(self):
        assert_equal(str(self.test_dna + self.test_dna), 'atgcatgc')
        assert_equal(str(self.test_dna.set_stranded('ss') +
                         self.test_dna.set_stranded('ss')), 'atgcatgc')

    def test_radd(self):
        assert_equal(str(sum([self.test_dna, self.test_dna])), 'atgcatgc')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_dna)

    def test_mul(self):
        assert_equal(str(self.test_dna * 4), 'atgcatgcatgcatgc')

        def mul_float(seq):
            return seq * 7.56

        assert_raises(TypeError, mul_float, self.test_dna)

        # TODO: reimplement this test using manual sequence input
        # def mul_incompatible(seq):
        #     return seq * 3

        # incompatible_seq = self.test_dna.copy()
        # incompatible_seq = incompatible_seq.five_resect(1)
        # incompatible_seq = incompatible_seq.reverse_complement()
        # incompatible_seq = incompatible_seq.five_resect(1)
        # incompatible_seq = incompatible_seq.reverse_complement()
        # assert_raises(Exception, mul_incompatible, incompatible_seq)

    def test_eq(self):
        assert_true(self.test_dna == sequence.DNA('atgc'))

    def test_ne(self):
        assert_true(self.test_dna != sequence.DNA('aagc'))

    def test_contains(self):
        assert_true('a' in self.test_dna)
        assert_true('t' in self.test_dna)
        assert_true('g' in self.test_dna)
        assert_true('c' in self.test_dna)
        assert_false('u' in self.test_dna)

    def test_flip(self):
        flipped = self.test_dna.flip()
        assert_equal(str(self.test_dna), flipped._bottom)
        assert_equal(self.test_dna._bottom, str(flipped))


def test_bad_bottom_init():
    def init_dna(top, bottom):
        sequence.DNA(top, bottom=bottom)

    assert_raises(ValueError, sequence.DNA, 'atgc', bottom='at')
    assert_raises(ValueError, sequence.DNA, 'atgc', bottom='gggg')
    try:
        sequence.DNA('atgc', bottom='gcat')
    except:
        assert False


def test_stranded_init():
    ss_dna = sequence.DNA('atgc', stranded='ss')
    assert_true(all([base == '-' for base in ss_dna._bottom]))

    ds_dna = sequence.DNA('atgc')
    assert_equal(str(ds_dna), sequence.utils.reverse_complement(ds_dna._bottom,
                                                                'dna'))


def test_stranded_complemented():
    ss_dna = sequence.DNA('atgc', stranded='ss')
    r_ss_dna = ss_dna.reverse_complement()
    assert_equal(r_ss_dna.top(), 'gcat')
    assert_equal(r_ss_dna.bottom(), '----')


class TestFeatures(object):
    '''Test features model using DNA object.'''
    def __init__(self):
        self.dna = sequence.DNA('ATGC') * 50
        self.apply_features()

    def apply_features(self):
        misc_feature = sequence.Feature('Misc Feature', 1, 20, 'misc_feature')
        misc_1_feature = sequence.Feature('Misc Feature', 1, 20,
                                          'misc_feature', strand=1)
        coding_feature = sequence.Feature('Coding Feature', 21, 40, 'CDS')
        primer_feature = sequence.Feature('Primer Feature', 41, 60, 'primer')
        promoter_feature = sequence.Feature('Promoter Feature', 61, 80,
                                            'promoter')
        terminator_feature = sequence.Feature('Terminator Feature', 81, 100,
                                              'terminator')
        rbs_feature = sequence.Feature('RBS Feature', 101, 120, 'RBS')
        origin_feature = sequence.Feature('Origin Feature', 121, 140,
                                          'rep_origin')
        utr3_feature = sequence.Feature("3'UTR Feature", 141, 160, "3'UTR")
        origin_feature2 = sequence.Feature('Origin Feature', 161, 180,
                                           'rep_origin')

        input_features = [misc_feature, misc_1_feature, coding_feature,
                          primer_feature, promoter_feature, terminator_feature,
                          rbs_feature, origin_feature, utr3_feature,
                          origin_feature2]
        self.dna = sequence.DNA(str(self.dna), features=input_features)

    def test_good_features(self):
        for feature in self.dna.features:
            assert_true(feature.copy() in self.dna.features)

    def test_bad_feature(self):
        assert_raises(ValueError, sequence.DNA, 'atgc', features='duck')
        assert_raises(ValueError, sequence.DNA, 'atgc', features=['duck'])
        assert_raises(ValueError, sequence.Feature, 'yEVenus', 0, 717, 'duck')

    def test_rev_comp(self):
        rev = self.dna.reverse_complement()
        for feature, rev_feature in zip(self.dna.features, rev.features):
            assert_not_equal(feature.strand, rev_feature.strand)
            assert_equal(len(self.dna) - feature.start, rev_feature.stop)
            assert_equal(len(self.dna) - feature.stop, rev_feature.start)

    def test_extract(self):
        test_utr3_feature = sequence.Feature('3\'UTR Feature', 0, 19, '3\'UTR')
        extracted = self.dna.extract('3\'UTR Feature')
        assert_equal(test_utr3_feature, extracted.features[0])
        assert_equal(str(extracted), 'tgcatgcatgcatgcatgc')

        assert_raises(ValueError, self.dna.extract, 'duck')
        assert_raises(ValueError, self.dna.extract, 'Origin Feature')

    def test_getitem(self):
        subsequence = self.dna[30:100]
        remaining_features = [sequence.Feature('Primer Feature', 11, 30,
                                               'primer'),
                              sequence.Feature('Promoter Feature', 31, 50,
                                               'promoter'),
                              sequence.Feature('Terminator Feature', 51, 70,
                                               'terminator')]

        assert_equal(subsequence.features, remaining_features)
        assert_false(self.dna[10].features)
        new_seq = sequence.DNA('ATGC',
                               features=[sequence.Feature('A', 0, 0,
                                                          'misc_feature')])
        assert_equal(new_seq[0].features[0],
                     sequence.Feature('A', 0, 0, 'misc_feature'))

    def test_delitem(self):
        copy = self.dna.copy()
        del copy[3]

        coding_feature = sequence.Feature('Coding Feature', 20, 39, 'CDS')
        primer_feature = sequence.Feature('Primer Feature', 40, 59, 'primer')
        promoter_feature = sequence.Feature('Promoter Feature', 60, 79,
                                            'promoter')
        terminator_feature = sequence.Feature('Terminator Feature', 80, 99,
                                              'terminator')
        rbs_feature = sequence.Feature('RBS Feature', 100, 119, 'RBS')
        origin_feature = sequence.Feature('Origin Feature', 120, 139,
                                          'rep_origin')
        utr3_feature = sequence.Feature('3\'UTR Feature', 140, 159, '3\'UTR')
        origin_feature2 = sequence.Feature('Origin Feature', 160, 179,
                                           'rep_origin')
        input_features_ref = [coding_feature, primer_feature,
                              promoter_feature, terminator_feature,
                              rbs_feature, origin_feature, utr3_feature,
                              origin_feature2]
        assert_equal(copy.features, input_features_ref)

    def test_ne(self):
        '''Test != operator'''
        assert_true(self.dna.features[0] != self.dna.features[4])
        assert_false(self.dna.features[0] != self.dna.features[0])


class TestRestrictionSite(object):
    '''Test RestrictionSite class.'''
    def __init__(self):
        self.ecorv = sequence.RestrictionSite(sequence.DNA('GATATC'), (3, 3),
                                              name='EcoRV')
        self.foki = sequence.RestrictionSite(sequence.DNA('GGATG'), (14, 18),
                                             name='FokI')

    def test_cuts_outside(self):
        '''Test cuts_outside method.'''
        assert_false(self.ecorv.cuts_outside())
        assert_true(self.foki.cuts_outside())

    def test_len(self):
        '''Test len function.'''
        assert_equal(len(self.ecorv), 6)
        assert_equal(len(self.foki), 5)
