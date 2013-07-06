'''
Test palindrome function from analysis module.

'''

from nose.tools import assert_true
from nose.tools import assert_false
from pymbt import analysis
from pymbt import sequence


def test_palindrome():
    palindromic_seq = sequence.DNA('ATGCGCAT')
    nonpalindromic_seq = sequence.DNA('ATGCGCAA')

    assert_true(analysis.palindrome(palindromic_seq))
    assert_false(analysis.palindrome(nonpalindromic_seq))
