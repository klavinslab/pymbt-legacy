import math
import time
from string import maketrans, translate
import multiprocessing
import csv

from pymbt.nupack import Nupack
from pymbt.oligo_synthesis.structure_windows import context_walk


def split_gene(seq,
               max_length=1100,
               core=60,
               context=90,
               step=10):

    eval_len = core + context
    # Trim down sequence to the necessary overlaps
    num_pieces = math.ceil(len(seq) / float(max_length))
    #TODO: 200 minimum window was chosen arbitrarily
    #       - how much wiggle room do we need?
    #      Should get a sense of how good/bad the score is, *then* decide
    #TODO: Should handle arbitrary number of overlaps, not just 2...
    if num_pieces * max_length - len(seq) < 200:
        num_pieces += 1
    if num_pieces == 1:
        seq_start_end = [(0, len(seq))]
    elif num_pieces == 2:
    #    seq_start_end = [(max_length,len(seq)-max_length)]
        seq_start_end = [(len(seq) - max_length, max_length)]
    elif num_pieces == 3:
        seq_start_end1 = (len(seq) - (2 * max_length - core), max_length)
        seq_start_end2 = (len(seq) - max_length, 2 * max_length - core)
        seq_start_end = [seq_start_end1, seq_start_end2]
    else:
        raise ValueError(
            'Sequences requiring greater than 3 pieces not yet implemented')

    seqs = [seq[x[0]:x[1]] for x in seq_start_end]
    walked = [context_walk(x, core, context, step, report=True) for x in seqs]

    scores = [[y[2] for y in x] for x in walked]
    maxscores = [x.index(max(x)) for x in scores]

    # TODO: fix redundant operations below - 'summary' is constructed twice
    max_summary = [walked[i][x] for i, x in enumerate(maxscores)]
    starts = [x[0] for x in seq_start_end]
    s_starts = [x[0] + starts[i] for i, x in enumerate(max_summary)]
    s_ends = [x[1] + starts[i] for i, x in enumerate(max_summary)]
    seq_scores = [x[2] for i, x in enumerate(max_summary)]
    seq_final = [seq[s_starts[i]:s_ends[i]] for i, x in enumerate(s_starts)]
    s1 = [seq[0:s_ends[i]] for i, x in enumerate(seqs)]
    s2 = [seq[s_starts[i]:] for i, x in enumerate(seqs)]

    #max_summary = [(seq_final[i], seq1[i], seq2[i], s_starts[i],
    #s_ends[i], seq_scores[i]) for i, x in enumerate(s_starts)]
    max_summary = []
    for i, x in enumerate(s_starts):
        a = (seq_final[i], s1[i], s2[i], s_starts[i], s_ends[i], seq_scores[i])
        max_summary.append(a)

    return(max_summary)

#TODO:
# want to ideally find global optimum configuration,
# but this becomes intractable w/ e.g. 100 scores+positions
# at 4 overlaps (5 pieces - future work).
# solution: cut it down to 50 right off the bat by throwing out half of the
# worst-scoring regions. Then, if we had 4 overlaps, we'd need
# to get 50^4 sums rather than 100^4, which is 16 times fewer
# (may need to employ more tricks than this as well
#   - trim by distance possibilities,e.g.)
# this trimming is a good idea - follows this general rule:
# 2 connections = sum_i_n a_i where a_i is just the ith integer e.g. 1,2,3,
# 3 connections = sum_i_n (a_i) + (a_i-1)


'''
# could compare to mfe
def run_mfe(sequence):
    mfe_np = Nupack(sequence,material='dna')
    mfe = mfe_np.mfe(0)
    mfe_np._close()
    return(mfe)

def mfe_score(seq):
    core_s = range(context-1,len(seq)-context-core,step)
    core_e = [ x+core for x in core_s ]
    cores = [seq[core_s[i]:core_e[i]] for i, x in enumerate(core_s)]
    total = len(cores)
    nupack_pool = multiprocessing.Pool()
    try:
        nupack_iterator = nupack_pool.imap(run_mfe,cores)
        # Watch progress
        while (True):
            completed = nupack_iterator._index
            if (completed == total):
                break
            else:
                print('(%s/%s) completed.') %(completed,total)
                time.sleep(4)
        mfes = [ x for x in nupack_iterator ]
        nupack_pool.close()
        nupack_pool.join()
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        nupack_pool.terminate()
        nupack_pool.close()
    return(mfes)
#mfes = mfe_score(seq)
#mfe_file = open('output_mfes.csv','wb')
#mfe_csv_writer = csv.writer(mfe_file,
                             delimiter=',',
                             quotechar="'",
                             quoting=csv.QUOTE_MINIMAL)
#mfe_csv_writer.writerow(['site','score'])
#for i in range(len(mfes)):
#    mfe_csv_writer.writerow([i+1,mfes[i]])
'''
