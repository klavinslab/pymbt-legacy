import time
import multiprocessing

from pymbt.dna_manipulation import reverse_complement
from pymbt.nupack import Nupack


def context_walk(seq, core_len, context_len, step, report=False):
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

    total = len(allseqs)
    nupack_pool = multiprocessing.Pool()
    msg = 'pair probabilities completed.'
    try:
        nupack_iterator = nupack_pool.imap(run_pairs, allseqs)
        # Watch progress
        t = 4
        while report:
            completed = nupack_iterator._index
            if (completed == total):
                break
            else:
                if t >= 4:
                    print('(%s/%s) ' % (completed, total) + msg)
                    t = 0
                t += 1
                time.sleep(1)
        all_pairs = [x for x in nupack_iterator]
        nupack_pool.close()
        nupack_pool.join()
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        nupack_pool.terminate()
        nupack_pool.close()

    allprobs = [x['probabilities'][-core_len:] for x in all_pairs]
    allscores = [sum(x) / len(x) for x in allprobs]

    # recondense the list
    lscores = allscores[0:len(allseqs) / 2]
    rscores = allscores[len(allseqs) / 2:]
    scores = [(lscores[i] + rscores[i]) / 2 for i, x in enumerate(lscores)]
    summary = []
    for i, x in enumerate(lseqs):
        summary.append((core_starts[i], core_ends[i], scores[i]))
        #print(seq[core_starts[i]:core_ends[i]])
        #print(lseqs[i])
        #print(rseqs[i])

    return(summary)


# Put nupack calculation in a function to enable
# parallel processing
def run_pairs(seq):
    np_run = Nupack(seq, 'dna')
    pairs = np_run.pairs(0)
    np_run._close()
    return(pairs)
