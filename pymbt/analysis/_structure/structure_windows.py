'''
Walk over a sequence in defined windows, evaluating for context-depenent
structure.

'''

from matplotlib import pylab
from pymbt.analysis import nupack_multiprocessing


class StructureWindows(object):
    '''
    Evaluate sequences in-context using Nupack's ppairs in windows along
    the sequence.

    '''

    def __init__(self, dna_object):
        '''
        :param dna_object: Input sequence - DNA object.
        :type seq: DNA object

        '''

        self.template = dna_object

        self.walked = []
        self.core_starts = []
        self.core_ends = []
        self.scores = []

    def run(self, window_size=60, context_len=90, step=10):
        '''
        Walk through the sequence of interest in windows of window_size,
        evaluate free (unbound) pair probabilities.

        :param window_size: Window size in base pairs.
        :type window_size: int
        :param context_len: The number of bases of context to use when
                            analyzing each window.
        :type context_len: int
        :param step: The number of base pairs to move for each new window.
        :type step: int

        '''

        self.walked = _context_walk(self.template, window_size, context_len,
                                    step)
        self.core_starts, self.core_ends, self.scores = zip(*self.walked)

        return self.walked

    def plot(self):
        '''
        Plot the results of ContextWalker.walk.

        '''
        if self.walked:
            fig = pylab.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(self.core_starts, self.scores, 'bo-')
            pylab.xlabel('Core sequence start position (base pairs).')
            pylab.ylabel('Score - Probability of being unbound.')
            pylab.show()
        else:
            raise Exception("Run calculate() first so there\'s data to plot!")


def _context_walk(dna_object, window_size, context_len, step):
    '''
    Generate context-dependent 'non-boundedness' series of scores for a given
    DNA sequence. Uses NUPACK's pair probabilities to derive the score.

    :param dna_object: DNA object.
    :type dna_object: DNA object
    :param window_size: Window size in base pairs.
    :type window_size: int
    :param context_len: The number of bases of context to use when analyzing
                        each window.
    :type context_len: int
    :param step: The number of base pairs to move for each new window.
    :type step: int

    '''

    # Generate window indices
    window_start_ceiling = len(dna_object) - context_len - window_size
    window_starts = range(context_len - 1, window_start_ceiling, step)
    window_ends = [start + window_size for start in window_starts]

    # Generate left-context subsequences
    l_starts = [step * i for i in range(len(window_starts))]
    l_seqs = []
    for start, end in zip(l_starts, window_ends):
        l_seqs.append(dna_object[start:end])

    # Generate right-context subsequences
    r_ends = [x + window_size + context_len for x in window_starts]
    r_seqs = []
    for start, end in zip(window_starts, r_ends):
        r_seqs.append(dna_object[start:end])
    r_seqs = [r_seq.reverse_complement() for r_seq in r_seqs]

    # Combine and calculate nupack pair probabilities
    seqs = l_seqs + r_seqs
    pairs_run = nupack_multiprocessing(seqs, 'dna', 'pairs', {'index': 0})
    # Focus on pair probabilities that matter - those in the window
    pairs = [run['probabilities'][-window_size:] for run in pairs_run]
    # Score by average pair probability
    lr_scores = [sum(pair) / len(pair) for pair in pairs]

    # Split into left-right contexts again and sum for each window
    l_scores = lr_scores[0:len(seqs) / 2]
    r_scores = lr_scores[len(seqs) / 2:]
    scores = [(l + r) / 2 for l, r in zip(l_scores, r_scores)]

    # Summarize and return window indices and score
    summary = []
    for start, end, score in zip(window_starts, window_ends, scores):
        summary.append((start, end, score))

    return summary