'''
Walk over a sequence in defined windows, evaluating for context-depenent
structure.

'''

from matplotlib import pylab
from pymbt.sequence_manipulation import reverse_complement
from pymbt.nupack import nupack_multiprocessing


class ContextWalker:
    '''
    Evaluate sequences in-context using Nupack's ppairs in windows along
    the sequence.

    '''

    def __init__(self, seq):
        '''
        :param seq: Input sequence.
        :type seq: str

        '''

        self.seq = seq
        self.walked = []
        self.core_starts = []
        self.core_ends = []
        self.scores = []

    def walk(self, core_len=60, context_len=90, step=10):
        '''
        Walk through the sequence of interest in windows of core_len, evaluate
        for free (unbound) pair probabilities.

        :param core_len: Window size in base pairs.
        :type core_len: int
        :param context_len: The number of bases of context to use when
                            analyzing each window.
        :type context_len: int
        :param step: The number of base pairs to move for each new window.
        :type step: int

        '''

        self.walked = context_walk(self.seq, core_len, context_len, step)
        self.core_starts = [x[0] for x in self.walked]
        self.core_ends = [x[1] for x in self.walked]
        self.scores = [x[2] for x in self.walked]

    def plot(self):
        '''Plot the results of ContextWalker.walk.'''
        if self.walked:
            fig = pylab.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(self.core_starts, self.scores, 'bo-')
            pylab.xlabel('Core sequence start position (base pairs).')
            pylab.ylabel('Score - Probability of being unbound.')
        else:
            print 'Run ContextWalker.calculate() first so there\'s data to \
                   \nplot!'


def context_walk(seq, core_len, context_len, step):
    '''
    Generate context-dependent 'non-boundedness' series of scores for a given
    DNA sequence. Uses NUPACK's pair probabilities to derive the score.

    :param seq: Input sequence.
    :type seq: str
    :param core_len: Window size in base pairs.
    :type core_len: int
    :param context_len: The number of bases of context to use when analyzing
                        each window.
    :type context_len: int
    :param step: The number of base pairs to move for each new window.
    :type step: int

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
    for i in range(len(lseqs)):
        summary.append((core_starts[i], core_ends[i], scores[i]))

    return summary
