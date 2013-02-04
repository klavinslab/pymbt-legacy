import time
import multiprocessing

from pymbt.dna_manipulation import reverse_complement
from pymbt.nupack import Nupack 

def context_walk(seq,core_len,context_len,step,report=False):
    # split into sections. not sophisticated
    # variable names are obscure
    core_starts = range(context_len-1,len(seq)-context_len-core_len,step)
    core_ends = [ x+core_len for x in core_starts ]
    lseq_starts = [ step*i for i,x in enumerate(core_starts) ]
    lseq_ends = core_ends
    rseq_starts = core_starts
    rseq_ends = [ x+core_len+context_len for x in rseq_starts ]
    lseqs = [ seq[lseq_starts[i]:lseq_ends[i]] for i,x in enumerate(lseq_starts) ]
    rseqs = [ seq[rseq_starts[i]:rseq_ends[i]] for i,x in enumerate(rseq_starts) ]
    rseqs = [ reverse_complement(x) for x in rseqs ]
    allseqs = lseqs + rseqs

    total = len(allseqs)
    nupack_pool = multiprocessing.Pool()
    try:
        nupack_iterator = nupack_pool.imap(run_pairs,allseqs)
        # Watch progress
        t = 4
        while (True) & report == True:
            completed = nupack_iterator._index 
            if (completed == total):
                break
            else:
                if t >= 4:
                    print('(%s/%s) completed.') %(completed,total)
                    t = 0
                t += 1
                time.sleep(1)
        all_pairs = [ x for x in nupack_iterator ]
        nupack_pool.close()
        nupack_pool.join()
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        nupack_pool.terminate()
        nupack_pool.close()

    allprobs = [ x['probabilities'][-core_len:] for x in all_pairs ]
    allscores = [ sum(x)/len(x) for x in allprobs ]

    # recondense the list
    lscores = allscores[0:len(allseqs)/2]
    rscores = allscores[len(allseqs)/2:]
    scores = [ (lscores[i]+rscores[i])/2 for i,x in enumerate(lscores) ]
    summary = [ (core_starts[i],core_ends[i],scores[i]) for i,x in enumerate(lseqs) ]

    return(summary)

# Put nupack calculation in a function to enable
# parallel processing
def run_pairs(seq):
    np_run = Nupack(seq,'dna')
    pairs = np_run.pairs(0)
    np_run._cleanup()
    return(pairs)
