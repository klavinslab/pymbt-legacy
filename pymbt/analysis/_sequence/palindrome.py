'''
Module for analyzing palindromic sequences.

'''


def palindrome(seq):
    '''
    Test whether a sequence is palindrome.

    :param seq: Sequence to analyze (DNA or RNA).
    :type seq: pymbt.sequence.DNA or pymbt.sequence.RNA

    '''

    seq_len = len(seq)
    wing = seq_len // 2

    if seq_len % 2 != 0:
        # Sequence has odd number of bases, need middle base to test complement
        l_wing = seq[0:wing + 1]
        r_wing = seq[wing:]
    else:
        # Sequence has even number of bases, can test non-overlapping seqs
        l_wing = seq[0: wing]
        r_wing = seq[wing:]

    if l_wing == r_wing.reverse_complement():
        return True
    else:
        return False
