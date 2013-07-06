'''
Module for analyzing palindromic sequences.

'''

import math
from pymbt import sequence


def palindrome(pattern):
    '''
    Check whether pattern is palindrome.
    :param pattern: pattern to test.
    :type pattern: str

    '''

    pattern = sequence.DNA(pattern)
    wing = int(math.floor(pattern / 2))

    if len(pattern) % 2 != 0:
        l_wing = pattern[0:wing + 1]
        r_wing = pattern[wing:]
    else:
        l_wing = pattern[0: wing]
        r_wing = pattern[wing:]
    if l_wing == r_wing.reverse_complement():
        return True
    else:
        return False
