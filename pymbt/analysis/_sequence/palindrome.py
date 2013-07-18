'''Analyze palindromic sequences.'''


def palindrome(seq):
    '''Test whether a sequence is palindrome.

    :param seq: Sequence to analyze (DNA or RNA).
    :type seq: pymbt.sequence.DNA or pymbt.sequence.RNA

    '''
    seq_len = len(seq)
    if seq_len % 2 == 0:
        # Sequence has even number of bases, can test non-overlapping seqs
        wing = seq_len / 2
        l_wing = seq[0: wing]
        r_wing = seq[wing:]
        if l_wing == r_wing.reverse_complement():
            return True
        else:
            return False
    else:
        # Sequence has odd number of bases and cannot be a palindrome
        return False
