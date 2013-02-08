import time
import random
import multiprocessing
import copy_reg
import types

import pymbt.orthoseq.OrthoSeq
from pymbt.nupack import Nupack


def monomers_concentration(sequence_list, mfe=True):
    # Run nupack's 'concentrations'
    nps = Nupack(sequence_list, 'dna')
    nc = nps.concentrations(2, conc=5e-7, mfe=mfe)
    concs = nc['concentration']
    types = nc['types']

    # Isolate the unbound monomer concentrations
    free_conc = sum([concs[i] for i, x in enumerate(types) if sum(x) == 1])
    free_fraction = free_conc / (2 * 5e-7)
    nps._close()  # Delete temp dir

    return(free_fraction)


def main(oligo_list, n=6, score_thresh=0.7):
    oligo_list_r = [OrthoSeq._revcomp(x) for x in oligo_list]
    mp = []
    for j in range(3):
        # Start with random oligo
        rc = random.choice(oligo_list)
        mp.append([oligo_list.pop(oligo_list.index(rc))])

        for i in range(1, n):
            np_list_ff = [mp[j] + [x] for x in oligo_list]
            np_list_fr = [mp[j] + [x] for x in oligo_list_r]
            np_list_both = np_list_ff + np_list_fr

            p = multiprocessing.Pool()
            monomers_iterator = p.imap(monomers_concentration, np_list_both)
            _progress(monomers_iterator, len(np_list_both), interval=5)
            monomers_list = [x for x in monomers_iterator]
            p.close()
            p.join()

            mon_f = [monomers_list.pop(0) for x in oligo_list]
            mon_r = monomers_list
            #TODO: expression below is checking min of the same thing?
            mon_mins = [min([x, x]) for i, x in enumerate(mon_f)]
            best_monomer_f = oligo_list.pop(mon_mins.index(max(mon_mins)))
            best_monomer_r = oligo_list_r.pop(mon_mins.index(max(mon_mins)))

            if i % 2 == 1:
                best_monomer = best_monomer_r
            else:
                best_monomer = best_monomer_f

            mp[j].append(best_monomer)

    # write a log
    f = open('last_run.txt', 'w')
    for j in mp:
        f.writelines('\n'.join(j))
        f.writelines('\n\n')
    f.close()
    return(mp)


def _progress(mp_iterator, total_jobs, interval=1):
    time_start = time.time()
    while (True):
        completed = mp_iterator._index
        if (completed == total_jobs):
            print '\n'
            break
        if completed > 0:
            time.sleep(interval)
            remaining = total_jobs - completed
            time_current = time.time()
            rate_current = float(completed) / (time_current - time_start)
            time_remaining = remaining / rate_current
            m = int(time_remaining) / 60  # minutes
            s = time_remaining - 60 * mins  # seconds
            print('Remaining: %(m).f min(s) %(s).f sec(s)' % {'m': m, 's': s})
            print('\r')
        else:
            time.sleep(.1)


# The following code enables pickling bound methods:
# Required for multiprocessing to operate
# on a method
def _pickle_method(method):
    name = method.__name__
    im_self = method.im_self
    im_class = method.im_class
    return(_unpickle_method, (name, im_self, im_class))


def _unpickle_method(func, im_self, im_class):
    return(getattr(im_self, func))

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
