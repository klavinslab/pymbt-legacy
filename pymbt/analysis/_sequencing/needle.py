'''Needleman-Wunsch alignment functions.'''
import multiprocessing
from pymbt import sequence
try:
    from calign import aligner, score_alignment
except ImportError:
    print "Failed to import cython aligner, so alignments will not work."
    pass


def needle(reference, query, gap_open=-15, gap_extend=0):
    '''Do a Needleman-Wunsch alignment.

    :param reference: Reference sequence.
    :type reference: pymbt.sequence.DNA
    :param query: Sequence to align against the reference.
    :type query: pymbt.sequence.DNA
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float
    :returns: (aligned reference, aligned query, score)
    :rtype: tuple of two pymbt.sequence.DNA instances and a float

    '''
    # Align using cython Needleman-Wunsch
    aligned_ref, aligned_res = aligner(str(reference),
                                       str(query),
                                       gap_open=gap_open,
                                       gap_extend=gap_extend,
                                       method="global_cfe",
                                       matrix="DNA_simple")

    # Score the alignment
    score = score_alignment(aligned_ref, aligned_res, gap_open, gap_extend,
                            "DNA_simple")

    return sequence.DNA(aligned_ref), sequence.DNA(aligned_res), score


def run_needle(args):
    """Run needle command using 4-tuple of the arguments (in the same order)
    as is used for needle. Necessary to make picklable function for
    multiprocessing."""
    return needle(*args)


def needle_multiprocessing(references, queries, gap_open=-15, gap_extend=0):
    """Batch process of sequencing split over several cores. Acts just like
    needle but sequence inputs are lists.

    :param references: References sequence.
    :type references: pymbt.sequence.DNA list
    :param queries: Sequences to align against the reference.
    :type queries: pymbt.sequence.DNA list
    :param gap_open: Penalty for opening a gap.
    :type gap_open: float
    :param gap_extend: Penalty for extending a gap.
    :type gap_extend: float
    :returns: a list of the same output as pymbt.sequence.needle
    :rtype: list

    """
    ## FIXME: remove this once issues figured out
    #return [needle(reference, query, gap_open=gap_open, gap_extend=gap_extend)
    #        for reference, query in zip(references, queries)]
    pool = multiprocessing.Pool()
    try:
        args_list = [list(x) + [gap_open, gap_extend] for x in
                     zip(references, queries)]
        aligned = pool.map(run_needle, args_list)
    except KeyboardInterrupt:
        pool.terminate()
        raise KeyboardInterrupt

    return aligned
