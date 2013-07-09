'''
Test palindrome function from analysis module.

'''

from nose.tools import assert_true
from nose.tools import assert_false
from pymbt import analysis
from pymbt import sequence


def test_palindrome():
    palindromic_seq_even = sequence.DNA('ATGCGCAT')
    nonpalindromic_seq_even = sequence.DNA('ATGCGCAA')

    assert_true(analysis.palindrome(palindromic_seq_even))
    assert_false(analysis.palindrome(nonpalindromic_seq_even))
