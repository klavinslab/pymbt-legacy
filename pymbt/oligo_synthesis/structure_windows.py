from pymbt.sequence_manipulation import reverse_complement
from pymbt.nupack import nupack_multiprocessing


def context_walk(seq, core_len, context_len, step, report=False):
    '''Generates context-dependent 'non-boundedness' series of scores for
       a given DNA sequence. Uses NUPACK's pair probabilities to derive
       the score.
           - seq:         input sequence.
           - core_len:    window size in base pairs.
           - context_len: the number of bases of context to use when analyzing
                          each window.
           - step:        the number of base pairs to move for each new window.
       '''
    # split into sections. not sophisticated
    # variable names are obscure
    adjusted = len(seq) - context_len - core_len
    core_starts = range(context_len - 1, adjusted, step)
    core_ends = [x + core_len for x in core_starts]
    l_starts = [step * i for i, x in enumerate(core_starts)]
    l_ends = core_ends
    r_starts = core_starts
    r_ends = [x + core_len + context_len for x in r_starts]
    lseqs = [seq[l_starts[i]:l_ends[i]] for i, x in enumerate(l_starts)]
    rseqs = [seq[r_starts[i]:r_ends[i]] for i, x in enumerate(r_starts)]
    rseqs = [reverse_complement(x) for x in rseqs]
    allseqs = lseqs + rseqs

    all_pairs = nupack_multiprocessing(allseqs, 'dna', 'pairs', {'strand': 0})
    allprobs = [x['probabilities'][-core_len:] for x in all_pairs]
    allscores = [sum(x) / len(x) for x in allprobs]

    # recondense the list
    lscores = allscores[0:len(allseqs) / 2]
    rscores = allscores[len(allseqs) / 2:]
    scores = [(lscores[i] + rscores[i]) / 2 for i, x in enumerate(lscores)]
    summary = []
    for i, x in enumerate(lseqs):
        summary.append((core_starts[i], core_ends[i], scores[i]))

    return summary
