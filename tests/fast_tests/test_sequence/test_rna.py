'''
Tests for the RNA sequence class.

'''

from pymbt import sequence
from nose.tools import assert_equal, assert_false, assert_true, assert_raises


class TestRNA(object):
    '''
    Testing class for sequence.RNA

    '''

    def __init__(self):
        self.test_rna = sequence.RNA('augc')

    def test_reverse_complement(self):
        assert_equal(self.test_rna.reverse_complement().top, 'gcau')

    def test_five_resect(self):
        assert_equal(self.test_rna.five_resect(2).top, 'gc')

    def test_three_resect(self):
        assert_equal(self.test_rna.three_resect(2).top, 'au')

    def test_locate(self):
        assert_equal(self.test_rna.locate('au'), [0])
        assert_equal(self.test_rna.locate('gc'), [2])
        assert_equal(len(self.test_rna.locate('augg')), 0)

    def test_copy(self):
        assert_equal(self.test_rna, self.test_rna.copy())

    def test_getitem(self):
        assert_equal(self.test_rna[0].top, 'a')
        assert_equal(self.test_rna[1].top, 'u')
        assert_equal(self.test_rna[2].top, 'g')
        assert_equal(self.test_rna[3].top, 'c')
        assert_equal(self.test_rna[-1].top, 'c')

    def test_delitem(self):
        copy0 = self.test_rna.copy()
        del copy0[0]
        assert_equal(copy0.top, 'ugc')
        copy1 = self.test_rna.copy()
        del copy1[1]
        assert_equal(copy1.top, 'agc')
        copy2 = self.test_rna.copy()
        del copy2[2]
        assert_equal(copy2.top, 'auc')
        copy3 = self.test_rna.copy()
        del copy3[3]
        assert_equal(copy3.top, 'aug')
        copy_1 = self.test_rna.copy()
        del copy_1[-1]
        assert_equal(copy_1.top, 'aug')

    def test_setitem(self):
        copy0 = self.test_rna.copy()
        copy0[0] = 'u'
        assert_equal(copy0.top, 'uugc')
        copy1 = self.test_rna.copy()
        copy1[1] = 'a'
        assert_equal(copy1.top, 'aagc')
        copy2 = self.test_rna.copy()
        copy2[2] = 'a'
        assert_equal(copy2.top, 'auac')
        copy3 = self.test_rna.copy()
        copy3[3] = 'a'
        assert_equal(copy3.top, 'auga')
        copy_1 = self.test_rna.copy()
        copy_1[-1] = 'a'
        assert_equal(copy_1.top, 'auga')

    def test_repr(self):
        expected_repr = 'RNA:\naugc\n----'
        assert_equal(repr(self.test_rna), expected_repr)

        repr_1 = 'RNA:\naugcaugcaugcaugcaugcaugcaugcaugcaugcaugc ... '
        repr_2 = 'augcaugcaugcaugcaugcaugcaugcaugcaugcaugc\n'
        repr_3 = '---------------------------------------- ... '
        repr_4 = '----------------------------------------'
        expected_long_repr = repr_1 + repr_2 + repr_3 + repr_4
        assert_equal(repr(self.test_rna * 50), expected_long_repr)

    def test_str(self):
        assert_equal(str(self.test_rna), 'augc')

    def test_len(self):
        assert_equal(len(self.test_rna), 4)

    def test_add(self):
        assert_equal((self.test_rna + self.test_rna).top, 'augcaugc')

    def test_radd(self):
        assert_equal(sum([self.test_rna, self.test_rna]).top, 'augcaugc')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_rna)

    def test_mul(self):
        assert_equal((self.test_rna * 4).top, 'augcaugcaugcaugc')

        def mul_float(seq):
            return seq * 7.56

        assert_raises(TypeError, mul_float, self.test_rna)

    def test_eq(self):
        assert_true(self.test_rna == sequence.RNA('augc'))

    def test_ne(self):
        assert_true(self.test_rna != sequence.RNA('aagc'))

    def test_contains(self):
        assert_true('a' in self.test_rna)
        assert_true('u' in self.test_rna)
        assert_true('g' in self.test_rna)
        assert_true('c' in self.test_rna)
        assert_false('t' in self.test_rna)
