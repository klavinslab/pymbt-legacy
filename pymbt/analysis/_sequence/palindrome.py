'''
Module for analyzing palindromic sequences.

'''


def palindrome(pattern):
    '''
    Check whether pattern is palindrome.
    :param pattern: pattern to test.
    :type pattern: DNA

    '''

    pattern_len = len(pattern)
    wing = pattern_len // 2

    if pattern_len % 2 != 0:
        l_wing = pattern[0:wing + 1]
        r_wing = pattern[wing:]
    else:
        l_wing = pattern[0: wing]
        r_wing = pattern[wing:]
    if l_wing == r_wing.reverse_complement():
        return True
    else:
        return False
