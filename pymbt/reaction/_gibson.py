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
        # validation

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
    # 1. Take the first sequence and find all matches
    # Get graphs:
    #   a) pattern watson : targets watson
    #   b) pattern watson : targets crick
    #   c) pattern crick: targets watson
    #   d) pattern crick: targets crick
    pattern = working_list[0]
    targets = working_list[1:]

    # Output graph nodes of terminal binders:
    #   (destination, size, strand1, strand2)
    def graph_strands(strand1, strand2):
        graph = []
        for i, target in enumerate(targets):
            matchlen = homology_report(pattern, target, strand1, strand2,
                                       cutoff=homology, min_tm=tm)
            if matchlen:
                graph.append((i, matchlen, strand1, strand2))
        return graph

    graph_ww = graph_strands("w", "w")
    graph_wc = graph_strands("w", "c")
    graph_cw = graph_strands("c", "w")
    graph_cc = graph_strands("c", "c")
    graphs_w = graph_ww + graph_wc
    graphs_c = graph_cw + graph_cc
    graphs = graphs_w + graphs_c

    # 2. See if there's more than one result on a strand.
    # If so, throw an exception.
    if len(graphs_w) > 1 or len(graphs_c) > 1:
        raise AmbiguousGibsonError('multiple compatible ends.')
    if len(graphs_w) == len(graphs_c) == 0:
        raise ValueError('Failed to find compatible Gibson ends.')

    # 3. There must be one result. Where is it?
    # If there's one result on each strand, go with the one that matches the
    # pattern watson strand (will occur first - index 0)
    match = graphs[0]

    # 4. Combine pieces together
    # 4a. Orient pattern sequence
    if match[2] == "c":
        left_side = pattern.reverse_complement()
    else:
        left_side = pattern
    # 4b. Orient target sequence
    if match[3] == "w":
        right_side = working_list.pop(match[0] + 1).reverse_complement()
    else:
        right_side = working_list.pop(match[0] + 1)

    working_list[0] = left_side + right_side[match[1]:]
    return working_list


def _fuse_last(working_list, homology, tm):
    '''With one sequence left, attempt to fuse it to itself.

    :param homology: length of terminal homology in bp.
    :type homology: int
    :raises: AmbiguousGibsonError if either of the termini are palindromic
             (would bind self-self).
             ValueError if the ends are not compatible.

    '''
    #1. Construct graph on self-self
    #   (destination, size, strand1, strand2)
    pattern = working_list[0]

    def graph_strands(strand1, strand2):
        matchlen = homology_report(pattern, pattern, strand1, strand2,
                                   cutoff=homology, min_tm=tm, top_two=True)
        if matchlen:
            # Ignore full-sequence matches
            # HACK: modified homology_report to accept top_two. It should
            # really just ignore full-length matches
            for length in matchlen:
                if length != len(pattern):
                    return (0, length, strand1, strand2)
        else:
            return []

    # cw is redundant with wc
    graph_ww = graph_strands("w", "w")
    graph_wc = graph_strands("w", "c")
    graph_cc = graph_strands("c", "c")
    if graph_ww + graph_cc:
        raise AmbiguousGibsonError("Self-self binding during circularization.")
    if not graph_wc:
        raise ValueError("Failed to find compatible ends for circularization.")

    working_list[0] = working_list[0][:-graph_wc[1]].circularize()

    return working_list


def homology_report(seq1, seq2, strand1, strand2, cutoff=0, min_tm=65.0,
                    top_two=False):
    '''Given two sequences (seq1 and seq2), report the size of all perfect
    matches between the 3' end of the top strand of seq1 and the 3' end of
    either strand of seq2. In short, in a Gibson reaction, what would bind the
    desired part of seq1, given a seq2?

    :param seq1: Sequence for which to test 3\' binding of a single strand to
                 seq2.
    :type seq1: pymbt.DNA
    :param seq2: Sequence for which to test 3\' binding of each strand to seq1.
    :type seq1: pymbt.DNA
    :param strand1: w (watson) or c (crick) - which strand of seq1 is being
                        tested.
    :type strand1ed: str
    :param strand2: w (watson) or c (crick) - which strand of seq2 is being
                        tested.
    :type strand2ed: str
    :param cutoff: size cutoff for the report - if a match is lower, it's
                   ignored
    :type cutoff: int
    :param min_tm: Minimum tm value cutoff - matches below are ignored.
    :type min_tm: float
    :returns: List of left and right identities.
    :rtype: list of ints

    '''
    # Go through each other sequence, both forward and revcomp sequences,
    # and check whether first base matches last / first of initial sequence
    if strand1 == "w":
        seq1_chunks = [seq1[-(i + 1):] for i in range(len(seq1))]
    else:
        # Inefficient (lots of reverse complementing)
        seq1_chunks = [seq1.reverse_complement()[-(i + 1):] for i in
                       range(len(seq1))]

    if strand2 == "w":
        seq2_chunks = [seq2[-(i + 1):] for i in range(len(seq2))]
    else:
        # Inefficient (lots of reverse complementing)
        seq2_chunks = [seq2.reverse_complement()[-(i + 1):] for i in
                       range(len(seq2))]

    target_matches = []
    for s1, s2 in zip(seq1_chunks, seq2_chunks):
        s1len = len(s1)
        # Inefficient! (reverse complementing a bunch of times)
        if s1.top() == s2.reverse_complement().top():
            if s1len >= cutoff:
                tm = pymbt.analysis.tm(s1)
                if tm >= min_tm:
                    target_matches.append(s1len)

    target_matches.sort()
    if not top_two:
        return 0 if not target_matches else target_matches[0]
    else:
        return 0 if not target_matches else target_matches[0:2]
