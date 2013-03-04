# Alternative version of OrthoSeq -
# Instead of evaluating individual pairwise interactions,
# look at the whole pool at once, every time - that way
# we get what we *really* need - sequences that can be
# placed into one big pot and not interact (or minimize that interaction)
#
# Generates a set of orthogonal DNA sequences for a given protein sequence
# Currently limited by the size of the protein sequence - ~10 residues max
# To increase size, need to optimize search method and take steps to reduce
# memory usage (write to files and iterate/queue)

#TODO: shouldn't pickle bound methods - use hidden functions (_function)
# necessary to pickle bound methods
# TODO: Implement signal handler so it exits cleanly on ctrl-c
#       Use 'try' approach
# TODO: Many parts of the while loop are repetitive
#       Roll into for loop or something.

import types
import copy_reg

import csv
from datetime import datetime
from itertools import combinations
from math import factorial
import multiprocessing
import os
import random
import socket
from string import maketrans, translate
import time

from pymbt.nupack import Nupack
from pymbt.sequence_generation import random_codons
from pymbt.sequence_generation import weighted_codons
from pymbt.sequence_manipulation import reverse_complement as r_c
from pymbt.sequence_manipulation import check_alphabet 
from pymbt.common_data import codon_freq_sc_nested

class OrthoSeq:
    def __init__(self,
                 prot_seq,
                 T=50,
                 min_free=0.95,
                 conc=5e-7,
                 nullinc_max=1e5):
        # Make sure sequence only includes amino acids
        check_alphabet(prot_seq, material='pep') 

        self.T = T
        self.prot_seq = prot_seq
        self.min_free = min_free
        self.conc = conc
        self.nullinc_max = nullinc_max

    def orthogonal_sequences(self,
                             n,
                             m=0,
                             resume=False,
                             wholemat=False,
                             stop=False,
                             report=True,
                             weighted=True,
                             freq_threshold=0.5):

        # TODO: 'stop' and 'report_nupack' may be deprecated
        '''n: is the final number of orthogonal sequences to output.
           m: is the number of new sequences to try out on each iteration
           resume: specifies a path to a dir from a previous run, and will cause
           orthogonal_sequences to pick up where it left off.
           wholemat: if True, m worst oligos are thrown out on each iteration
                     if False, single best is kept. There is a difference
                     but it should be better documented.
           stop: 
           report_nupack:
        '''
        self.weighted = weighted
        self.freq_threshold = freq_threshold
        # The default combination number should be done dynamically.
        # Solved combinations equation for m.
        # TODO: verify that this is true -  this is not m choose n
        n_combinations = factorial(n + m) / (factorial(m) * factorial(n))
        # Generate N oligos that meet min_free threshold
        if resume:
            oligo_list = resume
        else:
            i_pool = multiprocessing.Pool()
            min_f = self.min_free
            i_iter = i_pool.imap(self.gen_oligo, [min_f for x in range(n)])
            oligo_list = [x for x in i_iter]
            i_pool.close()
            i_pool.join()

        # Generate expected monomer concentrations matrices
        # for both forward-forward and forward-reverse
        c_mat_f = self.pairwise_monomer_concs(oligo_list)
        c_mat_r = self.pairwise_monomer_concs(oligo_list, reverse=True)

        # Initial loop takes too long for n > 3, so start at
        # 3 and add more as some are found?
        save_pos = []
        loop_count = 0
        time_elapsed = 0

        nullinc = 0
        last_min_free = 0
        while True:
            loop_start = time.time()
            loop_count += 1

            n_pool = multiprocessing.Pool()
            o_iterator = n_pool.imap(self.gen_oligo, [min_f for x in range(m)])
            new_oligos = [x for x in o_iterator]
            n_pool.close()
            n_pool.join()

            # Now calculate larger pool of oligos.
            # It is probably safe to only calculate:
            #   1. Each oligo in a 'forward' version against all
            #      other reversed oligos in a size-N pool
            #   2. Repeat this for each oligo
            #   3. Discard the m worst oligos, repeat

            oligo_list += new_oligos

            # Generate combinations - (m+n) choose n 
            combos = [list(x) for x in combinations(oligo_list, n)]

            # For each combination of length n,
            # make n lists of oligos where each one
            # is paired with the reverse of all the others
            # So each ith list contains every 'oligo vs reverse of all others'
            # possibility
            full_combos = [[] for i, x in enumerate(combos)]
            for i, v in enumerate(combos):
                for j in v:
                    newlist = [j]
                    newlist += [r_c(x) for x in v if x is not j]
                    full_combos[i].append(newlist)


            # Run that list through nupack
            o_pool = multiprocessing.Pool()
            o_pool_iter = o_pool.imap(self.n_p_nupack, full_combos)
            # Report progress
            if report:
                tot = len(full_combos)
                self._multiprocessing_progress(o_pool_iter, tot, interval=5)
            mon_concs = [x for x in o_pool_iter]
            o_pool.close()
            o_pool.join()

            # Flatten lists for easier scoring (min conc, mean conc)
            mon_concs_flat = [[w for v in x for w in v] for x in mon_concs]

            # Keep the best-scoring combo
            mc_mins = [min(x) for x in mon_concs_flat]
            best_combo = mc_mins.index(max(mc_mins))
            oligo_list = combos[best_combo]
            mon_concs_left = mon_concs_flat[best_combo]

            # Scoring - mean free conc, minimum free conc
            mean_free = (sum(mon_concs_left) / len(mon_concs_left)) / self.conc
            min_free = min(mon_concs_left) / self.conc
            best_concs = mon_concs[best_combo]
            met = [1 for x in best_concs if min(x) / self.conc >= self.min_free]
            met = sum(met)

            # Did the score improve? If not, increment counter
            if min_free == last_min_free:
                nullinc += 1
            else:
                nullinc = 0

            #######################################
            # Logging - provides data trail and allows resuming/passing along
            # an extended attempt
            #######################################

            # Initial file setup for logging
            if loop_count == 1:
                # Use date, time, host, and peptide sequence to uniquify file
                current_date = datetime.now().strftime('%Y%m%d--%H%M%S')
                host = socket.gethostname()
                cur_id = current_date + '-' + host + '-' + self.prot_seq
                cur_dir = os.getcwd()
                fileprefix = '%s/%s-oligo_calc_' % (cur_dir, cur_id)

                # Write out setup info file
                infofile = open(fileprefix + 'info.txt', 'w')
                info_csv = csv.writer(infofile, quoting=csv.QUOTE_MINIMAL)
                info_csv.writerow(['sequence', 'oligo_n', 'oligo_m'])
                info_csv.writerow([self.prot_seq, n, m])
                infofile.close()

                # Set up data log file
                datafile = open(fileprefix + 'data.txt', 'w')
                data_csv = csv.writer(datafile, quoting=csv.QUOTE_MINIMAL)
                data_csv.writerow(['loop', 'nullinc', 'mean_free', 'min_free', 'met', 'time'])
                datafile.close()

            # Write the current set of oligos to file (enables resuming)
            oligofile = open(fileprefix + 'latest.txt', 'w')
            for i in oligo_list:
                oligofile.write(i + '\n')
            oligofile.close()

            # Update the data log file
            datafile = open(fileprefix + 'data.txt', 'a')
            data_csv = csv.writer(datafile, quoting=csv.QUOTE_MINIMAL)
            time_diff = time.time() - loop_start
            time_elapsed += time_diff
            d = [loop_count, nullinc, mean_free, min_free, met, time_elapsed]
            data_csv.writerow(d)
            datafile.close()

            # Script stop decision:
            # Stop looping if either:
            #   1. the threshold is met
            #   2. the score has stopped improving by nullinc
            if met == len(oligo_list):
                break
            if nullinc == self.nullinc_max:
                break

            last_min_free = min_free

        return oligo_list

    def n_p_nupack(self, sequence_list):
        all_concs = []
        for x in sequence_list:
            n_np = Nupack(x, 'dna')
            n_concs = n_np.concentrations(2)['concentration']
            n_concs = n_concs[0:len(x)]
            all_concs.append(n_concs)
        return all_concs

    def gen_oligo(self, min_free):
        # Generate oligos until one has a self-self interaction below
        # the threshold.
        threshold_met = False 

        while not threshold_met:
            # Generate a randomly or weighted-randomly
            # coding sequence for the specified peptide
            # Ensure it has an mfe of 0
            if self.weighted:
                choice = weighted_codons.WeightedCodons(self.prot_seq,
                                                        frequency_table='sc',
                                                        material='pep')
                new_oligo = choice.generate()
            else:
                choice = random_codons.RandomCodons(self.prot_seq, threshold=self.freq_threshold)
                new_oligo = choice.generate()

            # Check oligo for self-self binding and mfe
            oligo_monomer = self.monomers_concentration(2 * [new_oligo])
            oligo_np = Nupack([new_oligo], 'dna')
            oligo_mfe = oligo_np.mfe(0)
            oligo_np.close()
            if oligo_monomer >= min_free and oligo_mfe == 0.0:
                threshold_met = True

        return new_oligo

    def monomers_concentration(self, sequence_list, mfe=True):
        # Run nupack's 'concentrations'
        np = Nupack(sequence_list, 'dna')
        nc = np.concentrations(2, conc=self.conc, mfe=mfe)
        concs = nc['concentration']
        types = nc['types']

        # Isolate the unbound monomer concentrations
        free_conc = sum([concs[i] for i, x in enumerate(types) if sum(x) == 1])
        free_fraction = free_conc / (2 * self.conc)
        np.close()  # Delete temp dir

        return free_fraction

    def pairwise_monomer_concs(self, seq_list, reverse=False):
        # Generate upper triangular matrix of sequence pairs.

        # If calculating forward-reverse concentrations, ignore
        # self-self as it will always interact strongly

        sl = seq_list 

        seq_pairs = []
        if reverse:
            for i, x in enumerate(sl):
                for j in range(i, len(sl)):
                    if i is not j:
                        seq_pairs.append([sl[i], r_c(sl[j])])
        else:
            for i, x in enumerate(sl):
                for j in range(i, len(sl)):
                    seq_pairs.append([sl[i], sl[j]])

        # Calculate unbound monomer concentrations for all pairs
        p = multiprocessing.Pool()
        pairwise_iterator = p.imap(self.monomers_concentration, seq_pairs)
        p_list = [x for x in pairwise_iterator]
        p.close()
        p.join()

        # Convert results (1-dimensional list) into upper triangular matrix
        # (list of lists)
        cm = []  # concentrations matrix
        for i, x in enumerate(sl):
            if reverse:
                # [] is used as placeholder for self-self when calculating
                # forward-reverse unbound monomer concentrations
                cur_bl = [[] for x in range(i + 1)]
                iter_range = range((len(sl) - i - 1))
                cur_row = cur_bl + [p_list.pop(0) for j in iter_range]
            else:
                cur_bl = [[] for x in range(i)]
                iter_range = range((len(sl) - i))
                cur_row = cur_bl + [p_list.pop(0) for j in iter_range]
            cm.append(cur_row)

        # Fill in the rest of the matrix
        for i, x in enumerate(cm):
            for j, y in enumerate(x):
                cm[j][i] = y

        return cm

    def _multiprocessing_progress(self, mp_iterator, total_jobs, interval=1):
        time_start = time.time()
        while (True):
            completed = mp_iterator._index
            if (completed == total_jobs):
                print '\n'
                break
            if completed > 0 & completed % 20 == 0:
                remaining = total_jobs - completed
                time_current = time.time()
                rate_current = float(completed) / (time_current - time_start)
                time_remaining = remaining / rate_current
                m = int(time_remaining) / 60  # minutes
                s = time_remaining - 60 * m  # seconds
                print("Remaining: %s:%.2f min(s)\r" % (m, s))
                time.sleep(interval)
            else:
                time.sleep(.1)

# The following code enables pickling bound methods:
# Required for multiprocessing to function
def _pickle_method(method):
    name = method.__name__
    im_self = method.im_self
    im_class = method.im_class
    return _unpickle_method, (name, im_self, im_class)


def _unpickle_method(func, im_self, im_class):
    return getattr(im_self, func)


copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


if __name__ == "__main__":
    import ConfigParser

    # Read config 
    current_path = os.path.abspath(os.path.dirname(__file__))
    config_path = current_path + '/orthogonal_sequences.cfg'
    config = ConfigParser.RawConfigParser()
    config.read(config_path)

    # Parse config into parameters
    sequence = config.get('Main Settings', 'protein sequence')
    n = config.getint('Main Settings', 'n sequences')
    m = config.getint('Main Settings', 'm sequences')
    T = config.getfloat('Main Settings', 'temperature')
    min_free = config.getfloat('Main Settings', 'min unbound monomer score')
    con = config.getfloat('Main Settings', 'species concentration')
    wholemat = config.getboolean('Main Settings', 'whole matrix')
    nm = config.getfloat('Main Settings', 'no increase limit')
    resume = config.get('Main Settings', 'resume')

    # Run OrthoSeq.orthogonal_sequences
    OS = OrthoSeq(sequence, T=T, min_free=min_free, conc=con, nullinc_max=nm)
    if resume == 'False':
        OS.orthogonal_sequences(n, m=m, wholemat=wholemat)
    else:
        resfile = open(resume, 'r')
        resume = resfile.readlines()
        resfile.close()
        resume = [x.strip() for x in resume]
        OS.orthogonal_sequences(n, m=m, wholemat=wholemat, resume=resume)
