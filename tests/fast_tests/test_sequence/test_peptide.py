'''
Tests for the Peptide sequence class.

'''

from pymbt import sequence
from nose.tools import assert_equal, assert_false, assert_true, assert_raises


class TestPeptide(object):
    '''
    Testing class for sequence.Peptide

    '''

    def __init__(self):
        self.test_peptide = sequence.Peptide('mkgp')

    def test_locate(self):
        assert_equal(self.test_peptide.locate('mk'), [0])
        assert_equal(self.test_peptide.locate('gp'), [2])
        assert_equal(len(self.test_peptide.locate('augg')), 0)

    def test_copy(self):
        assert_equal(self.test_peptide, self.test_peptide.copy())

    def test_getitem(self):
        assert_equal(str(self.test_peptide[0]), 'm')
        assert_equal(str(self.test_peptide[1]), 'k')
        assert_equal(str(self.test_peptide[2]), 'g')
        assert_equal(str(self.test_peptide[3]), 'p')
        assert_equal(str(self.test_peptide[-1]), 'p')

    def test_delitem(self):
        copy0 = self.test_peptide.copy()
        del copy0[0]
        assert_equal(str(copy0), 'kgp')
        copy1 = self.test_peptide.copy()
        del copy1[1]
        assert_equal(str(copy1), 'mgp')
        copy2 = self.test_peptide.copy()
        del copy2[2]
        assert_equal(str(copy2), 'mkp')
        copy3 = self.test_peptide.copy()
        del copy3[3]
        assert_equal(str(copy3), 'mkg')
        copy_1 = self.test_peptide.copy()
        del copy_1[-1]
        assert_equal(str(copy_1), 'mkg')

    def test_setitem(self):
        copy0 = self.test_peptide.copy()
        copy0[0] = 'q'
        assert_equal(str(copy0), 'qkgp')
        copy1 = self.test_peptide.copy()
        copy1[1] = 'q'
        assert_equal(str(copy1), 'mqgp')
        copy2 = self.test_peptide.copy()
        copy2[2] = 'q'
        assert_equal(str(copy2), 'mkqp')
        copy3 = self.test_peptide.copy()
        copy3[3] = 'q'
        assert_equal(str(copy3), 'mkgq')
        copy_1 = self.test_peptide.copy()
        copy_1[-1] = 'q'
        assert_equal(str(copy_1), 'mkgq')

    def test_str(self):
        assert_equal(str(self.test_peptide), 'mkgp')

    def test_len(self):
        assert_equal(len(self.test_peptide), 4)

    def test_add(self):
        assert_equal(str((self.test_peptide + self.test_peptide)),
                     'mkgpmkgp')

    def test_radd(self):
        assert_equal(str(sum([self.test_peptide, self.test_peptide])),
                     'mkgpmkgp')

        def radd_800(seq):
            return 800 + seq

        assert_raises(TypeError, radd_800, self.test_peptide)

    def test_mul(self):
        assert_equal(str(self.test_peptide * 4), 'mkgpmkgpmkgpmkgp')

        def mul_float(seq):
            return seq * 7.56

        assert_raises(TypeError, mul_float, self.test_peptide)

    def test_eq(self):
        assert_true(self.test_peptide == sequence.Peptide('mkgp'))

    def test_ne(self):
        assert_true(self.test_peptide != sequence.Peptide('mkqp'))

    def test_contains(self):
        assert_true('m' in self.test_peptide)
        assert_true('k' in self.test_peptide)
        assert_true('g' in self.test_peptide)
        assert_true('p' in self.test_peptide)
        assert_false('q' in self.test_peptide)
