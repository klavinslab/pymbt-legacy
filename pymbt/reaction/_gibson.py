'''Gibson reaction simulation.'''
import pymbt.analysis


class AmbiguousGibsonError(ValueError):
    '''Exception to raise when Gibson is ambiguous.'''
    pass


def gibson(seq_list, linear=False, homology=10, tm=65.0):
        '''Simulate a Gibson reaction.

        :param seq_list: list of DNA sequences to Gibson
        :type seq_list: list of pymbt.DNA
        :param linear: Attempt to produce linear, rather than circular,
                       fragment from input fragments.
        :type linear: bool
        :param homology_min: minimum bp of homology allowed
        :type homology_min: int
        :param tm: Minimum tm of overlaps
        :type tm: float
        :returns: pymbt.reaction.Gibson instance.
        :raises: ValueError if any input sequences are circular DNA.

        '''
        # FIXME: Preserve features in overlap
        # TODO: set a max length?
        # TODO: add 'expected' keyword argument somewhere to automate

        # FIXME: why?
        # Remove any redundant (identical) sequences
        seq_list = list(set(seq_list))

        for seq in seq_list:
            if seq.topology == "circular":
                raise ValueError("Input sequences must be linear.")

        # Copy input list
        working_list = seq_list[:]
        # Attempt to fuse fragments together until only one is left
        while len(working_list) > 1:
            working_list = _find_fuse_next(working_list, homology, tm)
        if not linear:
            # Fuse the final fragment to itself
            working_list = _fuse_last(working_list, homology, tm)
        return working_list[0]


def _find_fuse_next(working_list, homology, tm):
    '''Find the next sequence to fuse, and fuse it (or raise exception).

    :param homology: length of terminal homology in bp
    :type homology: int
    :raises: AmbiguousGibsonError if there is more than one way for the
             fragment ends to combine.

    '''
    # 1. Analyze all non-first sequences for matches
    pattern = working_list[0]
    targets = working_list[1:]
    report = [homology_report(pattern, target) for target in targets]
    print report

    # 2. Throw out matches that are too short or have too low of tm
    for i, target in enumerate(report):
        for j, side in enumerate(target):
            # If too short, throw out match
            report[i][j] = [x for x in side if x > homology]
            # If tm is too low, throw out match
            if j == 0:
                seqs = [targets[i][:x] for x in side]
                report[i][j] = [x for k, x in enumerate(side) if
                                pymbt.analysis.tm(seqs[k]) > tm]
            else:
                seqs = [targets[i][-x:] for x in side]
                report[i][j] = [x for k, x in enumerate(side) if
                                pymbt.analysis.tm(seqs[k]) > tm]

    # 3. See if there's more than one result. If so, throw exception
    count = len([z for x in report for y in x for z in y])
    if not count:
        raise ValueError('Failed to find compatible Gibson ends.')
    elif count > 1:
        raise AmbiguousGibsonError('multiple compatible ends.')

    # 4. There must be one result. Where is it?
    for i, target in enumerate(report):
        for j, side in enumerate(target):
            if side:
                match = (i, j, side[0])

    # 5. Combine pieces together
    left_side = working_list[0]
    right_side = working_list.pop(match[0] + 1)
    if match[1] == 1:
        right_side = right_side.reverse_complement()

    working_list[0] = left_side + right_side[match[2]:]
    return working_list


def _fuse_last(working_list, homology, tm):
    '''With one sequence left, attempt to fuse it to itself.

    :param homology: length of terminal homology in bp.
    :type homology: int
    :raises: AmbiguousGibsonError if either of the termini are palindromic
             (would bind self-self).
             ValueError if the ends are not compatible.

    '''
    # 1. Get report on self-self
    pattern = working_list[0]
    report = homology_report(pattern, pattern)
    report_l = [x for x in report[0] if x != len(pattern)]
    report_l = [x for x in report_l if x > homology]
    report_l = [x for x in report_l if
                pymbt.analysis.tm(pattern[:x + 1]) > tm]
    if not report_l:
        raise ValueError('Failed to find compatible Gibson ends.')
    elif len(report_l) > 1:
        raise AmbiguousGibsonError('multiple compatible ends.')
    working_list[0] = pattern[:-report_l[0]].circularize()

    return working_list


def homology_report(seq1, seq2):
    '''Given two sequence inputs, find 3\' revcomp identities.

    :param seq1: Sequence to test 3\' end of Watson strand of
    :type seq1: pymbt.DNA
    :param seq2: Sequence to test both strands of, 3\' end
    :type seq1: pymbt.DNA
    :returns: List of left and right identities.
    :rtype: list of ints

    '''
    # Go through each other sequence, both forward and revcomp sequences,
    # and check whether first base matches last / first of initial sequence
    reference = str(seq1)
    watson_3 = str(seq2)
    crick_3 = str(seq2.reverse_complement())

    def right_matches_left(seq1, seq2, n):
        return seq1[-(n + 1):] == seq2[:n + 1]

    left_list = []
    for i in range(min(len(reference), len(watson_3))):
        if right_matches_left(reference, watson_3, i):
            left_list.append((i + 1))

    right_list = []
    for i in range(min(len(reference), len(crick_3))):
        if right_matches_left(reference, crick_3, i):
            right_list.append((i + 1))

    return [left_list, right_list]
