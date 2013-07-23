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
        assert_equal(self.test_dna.reverse_complement().top, 'gcat')

    def test_linearize(self):
        circ_dna = self.test_dna.circularize()
        assert_equal(circ_dna.linearize().top, 'atgc')
        assert_equal(circ_dna.linearize(0).top, 'atgc')
        assert_equal(circ_dna.linearize(1).top, 'tgca')
        assert_equal(circ_dna.linearize(2).top, 'gcat')
        assert_equal(circ_dna.linearize(3).top, 'catg')
        assert_equal(circ_dna.linearize(-1).top, 'catg')

        def linearize_linear(seq):
            return seq.linearize()
        assert_raises(ValueError, linearize_linear, circ_dna.linearize())

    def test_set_stranded(self):
        assert_equal(self.test_dna.set_stranded('ds'), self.test_dna)
        ss_dna = self.test_dna.set_stranded('ss')
        assert_equal(ss_dna.stranded, 'ss')
        ds_dna = self.test_dna.set_stranded('ds')
        assert_equal(ds_dna.stranded, 'ds')
        assert_equal(sequence.utils.reverse_complement(ds_dna.bottom, 'dna'),
                     ss_dna.top)

        r_ss_dna = ds_dna  # TODO: make sure this is worth testing
        r_ds_dna = self.test_dna.set_stranded('ds')
        assert_equal(r_ds_dna.reverse_complement().top, r_ss_dna.bottom)

        ds_to_ss_to_ds = self.test_dna.set_stranded('ss').set_stranded('ds')
        assert_equal(self.test_dna, ds_to_ss_to_ds)

        empty_top = reaction.three_resect(self.test_dna, 400)
        assert_equal(empty_top.set_stranded('ds'), self.test_dna)

        def bad_argument(seq):
            seq.set_stranded('duck')
        assert_raises(ValueError, bad_argument, self.test_dna)

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
        assert_equal(self.test_dna[0].top, 'a')
        assert_equal(self.test_dna[1].top, 't')
        assert_equal(self.test_dna[2].top, 'g')
        assert_equal(self.test_dna[3].top, 'c')
        assert_equal(self.test_dna[-1].top, 'c')

    def test_delitem(self):
        copy0 = self.test_dna.copy()
        del copy0[0]
        assert_equal(copy0.top, 'tgc')
        copy1 = self.test_dna.copy()
        del copy1[1]
        assert_equal(copy1.top, 'agc')
        copy2 = self.test_dna.copy()
        del copy2[2]
        assert_equal(copy2.top, 'atc')
        copy3 = self.test_dna.copy()
        del copy3[3]
        assert_equal(copy3.top, 'atg')
        copy_1 = self.test_dna.copy()
        del copy_1[-1]
        assert_equal(copy_1.top, 'atg')

    def test_setitem(self):
        copy0 = self.test_dna.copy()
        copy0[0] = 't'
        assert_equal(copy0.top, 'ttgc')
        copy1 = self.test_dna.copy()
        copy1[1] = 'a'
        assert_equal(copy1.top, 'aagc')
        copy2 = self.test_dna.copy()
        copy2[2] = 'a'
        assert_equal(copy2.top, 'atac')
        copy3 = self.test_dna.copy()
        copy3[3] = 'a'
        assert_equal(copy3.top, 'atga')
        copy_1 = self.test_dna.copy()
        copy_1[-1] = 'a'
        assert_equal(copy_1.top, 'atga')

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
        assert_equal((self.test_dna + self.test_dna).top, 'atgcatgc')
        assert_equal((self.test_dna.set_stranded('ss') +
                      self.test_dna.set_stranded('ss')).top, 'atgcatgc')

    def test_radd(self):
        assert_equal(sum([self.test_dna, self.test_dna]).top, 'atgcatgc')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_dna)

    def test_mul(self):
        assert_equal((self.test_dna * 4).top, 'atgcatgcatgcatgc')

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


def test_bad_bottom_init():
    def init_dna(top, bottom):
        sequence.DNA(top, bottom=bottom)

    assert_raises(ValueError, init_dna, 'atgc', 'at')
    assert_raises(ValueError, init_dna, 'atgc', 'gggg')
    try:
        sequence.DNA('atgc', bottom='gcat')
    except:
        assert False


def test_stranded_init():
    ss_dna = sequence.DNA('atgc', stranded='ss')
    assert_true(all([base == '-' for base in ss_dna.bottom]))

    ds_dna = sequence.DNA('atgc')
    assert_equal(ds_dna.top, sequence.utils.reverse_complement(ds_dna.bottom,
                                                               'dna'))


def test_stranded_complemented():
    ss_dna = sequence.DNA('atgc', stranded='ss')
    r_ss_dna = ss_dna.reverse_complement()
    assert_equal(ss_dna.bottom, sequence.utils.reverse_complement(r_ss_dna.top,
                                                                  'dna'))
    assert_equal(ss_dna.top, sequence.utils.reverse_complement(r_ss_dna.bottom,
                                                               'dna'))


def test_feature():
    def badtype():
        sequence.Feature('yEVenus', 0, 717, 'duck')
    assert_raises(ValueError, badtype)


class TestFeatures(object):
    '''Test features model using DNA object.'''
    def __init__(self):
        self.dna = sequence.DNA('ATGC') * 50
        self.apply_features()

    def apply_features(self):
        misc_feature = sequence.Feature('Misc Feature', 1, 20, 'misc')
        misc_1_feature = sequence.Feature('Misc Feature', 1, 20, 'misc',
                                          strand=1)
        coding_feature = sequence.Feature('Coding Feature', 21, 40, 'coding')
        primer_feature = sequence.Feature('Primer Feature', 41, 60, 'coding')
        promoter_feature = sequence.Feature('Promoter Feature', 61, 80,
                                            'promoter')
        terminator_feature = sequence.Feature('Terminator Feature', 81, 100,
                                              'terminator')
        rbs_feature = sequence.Feature('RBS Feature', 101, 120, 'rbs')
        origin_feature = sequence.Feature('Origin Feature', 121, 140, 'origin')
        utr3_feature = sequence.Feature('3\'UTR Feature', 141, 160, '3\'utr')

        input_features = [misc_feature, misc_1_feature, coding_feature,
                          primer_feature, promoter_feature, terminator_feature,
                          rbs_feature, origin_feature, utr3_feature]
        self.dna = sequence.DNA(self.dna.top, features=input_features)

    def test_good_features(self):
        for feature in self.dna.features:
            assert_true(feature.copy() in self.dna.features)

    def test_bad_feature(self):
        def duckfeature():
            sequence.DNA('atgc', features='duck')
        assert_raises(ValueError, duckfeature)

        def duckfeaturelist():
            sequence.DNA('atgc', features=['duck'])
        assert_raises(ValueError, duckfeaturelist)

    def test_rev_comp(self):
        rev = self.dna.reverse_complement()
        for feature, rev_feature in zip(self.dna.features, rev.features):
            assert_not_equal(feature.strand, rev_feature.strand)
            assert_equal(len(self.dna) - feature.start, rev_feature.stop)
            assert_equal(len(self.dna) - feature.stop, rev_feature.start)
