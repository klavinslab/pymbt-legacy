'''Generates a set of orthogonal DNA sequences for a given protein sequence
   Currently limited by the size of the protein sequence - ~10 residues max
   To increase size, need to optimize search method and take steps to reduce
   memory usage (write to files and iterate/queue)'''

# TODO: Many parts of the while loop are repetitive

import csv
import collections
from datetime import datetime
from itertools import combinations
import os
import socket
import time

from pymbt.nupack import nupack_multiprocessing
from pymbt.nupack import Nupack
from pymbt.sequence_generation import random_codons
from pymbt.sequence_generation import weighted_codons
from pymbt.sequence_manipulation import reverse_complement as r_c
from pymbt.sequence_manipulation import check_alphabet


class OrthoSeq:
    def __init__(self,
                 prot_seq,
                 T=50,
                 min_score=0.95,
                 n_attempt=1e5):
        # Make sure sequence only includes amino acids
        check_alphabet(prot_seq, material='pep')
        # Input sequence - a peptide
        self.prot_seq = prot_seq
        # Temp at which to evaluate all reactions
        self.T = T
        # Primary scoring metric - minimum average free monomer concentrations
        # - 0 would be fully-bound, 1 would be nothing bound
        self.min_score = min_score
        # Secondary scoring metric - stops the simulation after n
        # attempts that fail to increase the primary score
        self.n_attempt = n_attempt

    def start(self,
              n,
              m=0,
              resume=False,
              wholemat=False,
              report=True,
              weighted=True,
              freq_threshold=0.5):

        '''n: is the final number of orthogonal sequences to output.
           m: is the number of new sequences to try out on each iteration.
           resume: specifies a path to a dir from a previous run, and will
           cause start() to pick up where it left off.
           wholemat: if True, m worst oligos are thrown out on each iteration.
                     if False, single best is kept. There is a difference
                     but it should be better documented.
           weighted: boolean for whether codons should be weighted.
           freq_threshold: relative frequency threshold below which codons
                     should not be used.
        '''
        # Concentrations at which simulations are done
        conc = 5e-7
        # Arguments for _gen_oligo sequence generation
        oligo_args = (self.prot_seq,
                      self.min_score,
                      freq_threshold,
                      weighted,
                      conc)
        # Generate n oligos that meet score threshold
        # If there's a list of oligos to resume improving, start using it
        if resume:
            if len(resume) != n:
                raise Exception('Number of sequences to resume doesn\'t match\
                                 setting of \'n\'.')
            for x in resume:
                if type(x) != str:
                    raise Exception('\'resume\' input isn\'t strings!')
            oligo_list = resume
        # If not, generate initial set of oligos
        else:
            oligo_list = [_gen_oligo(*oligo_args) for x in range(n)]
            oligo_list = remove_redundant(oligo_list, oligo_args)
        time_start = time.time()
        loop_count = 0
        n_attempted = 0
        last_min_score = 0
        # central loop of the program - keeps running until either the
        # primary score is met or the score fails to increase for n attempts
        while True:
            loop_count += 1
            # Generate new oligo(s) to test in combination with the rest
            new_oligos = [_gen_oligo(*oligo_args) for x in range(m)]
            oligo_list += new_oligos
            # If two or more oligos are the same, replace with unique ones
            oligo_list = remove_redundant(oligo_list, oligo_args)
            # Generate all combinations - (m+n) choose n
            combos = [list(x) for x in combinations(oligo_list, n)]

            # For every combination of length m, generate m new lists where
            # all but one of the oligos is reverse-complemented.
            # reversed_combos is a list of lists of sequences, each entry
            # compatible with the Nupack module's complexes method
            reversed_combos = []
            for combo in combos:
                for seq in combo:
                    newlist = [seq]
                    newlist += [r_c(x) for x in combo if x != seq]
                    reversed_combos.append(newlist)

            # a list of Nupack's 'complexes' that's m * (m + n) long
            np_concs_raw = nupack_multiprocessing(reversed_combos,
                                                  'dna',
                                                  'concentrations',
                                                  {'max_complexes': 2})
            # just the concentrations list from each entry
            np_concs = [x['concentration'] for x in np_concs_raw]
            # just the types list from each entry
            np_types = [x['types'] for x in np_concs_raw]
            # the monomer concentrations from each entry
            mon_concs = []
            for i, x in enumerate(np_types):
                # if sum of the list is 1, that means it's a monomer conc
                # (e.g. [1, 0, 0, 0])
                temp_list = []
                for j, v in enumerate(x):
                    if sum(v) == 1:
                        temp_list.append(np_concs[i][j])
                mon_concs.append(temp_list)

            # Group by combo
            grouped_mon_concs = unflatten(mon_concs, m)

            # Find the worst binder within each of the m
            # rev-comp variations per combo
            grouped_mins = [[min(y) for y in x] for x in grouped_mon_concs]

            # Find the worst binder for each combination
            mon_conc_mins = [min(x) for x in grouped_mins]
            best_index = mon_conc_mins.index(max(mon_conc_mins))
            best_concs = grouped_mins[best_index]
            # Keep the best-scoring combination of sequences
            best_combo = combos[best_index]

            # Scoring - mean min free conc, minimum free conc
            mean_score = (sum(best_concs) / len(best_concs)) / conc
            min_score = min(best_concs) / conc
            met = 0
            for x in best_concs:
                if (x / conc) >= self.min_score:
                    met += 1

            # Did the minimum score improve? If not, increment counter
            if min_score == last_min_score:
                n_attempted += 1
            else:
                n_attempted = 0
                oligo_list = best_combo

            # Record how much time has passed
            timer = time.time() - time_start

            # Log
            self.log(loop_count,
                     n_attempted,
                     mean_score,
                     min_score,
                     met,
                     timer,
                     oligo_list)

            # Script stop decision:
            # Stop looping if either:
            #   1. the threshold is met
            #   2. the minimum score has stopped improving by n_attempt
            if met == len(oligo_list):
                break
            if n_attempted == self.n_attempt:
                break

            last_min_score = min_score

        return oligo_list

    def log(self,
            loop_count,
            n_attempted,
            mean_score,
            min_score,
            met,
            timer,
            oligo_list):
        '''provides data trail and allows resuming/passing along
           an extended attempt'''

        # Initial file setup for logging
        if loop_count == 1:
            # Use date, time, host, and peptide sequence to uniquify file
            current_date = datetime.now().strftime('%Y%m%d--%H%M%S')
            host = socket.gethostname()
            cur_id = current_date + '-' + host + '-' + self.prot_seq
            cur_dir = os.getcwd()
            file_prefix = '%s/%s-oligo_calc_' % (cur_dir, cur_id)

            # Write out setup info file
            info_file = open(file_prefix + 'info.txt', 'w')
            info_csv = csv.writer(info_file, quoting=csv.QUOTE_MINIMAL)
            info_csv.writerow(['sequence', 'oligo_n', 'oligo_m'])
            info_csv.writerow([self.prot_seq, n, m])
            info_file.close()

            # Set up data log file
            data_file = open(file_prefix + 'data.csv', 'w')
            data_csv = csv.writer(data_file, quoting=csv.QUOTE_MINIMAL)
            cols = ['loop',
                    'n_attempted',
                    'mean_score',
                    'min_score',
                    'met',
                    'time']
            data_csv.writerow(cols)
            data_file.close()

        # Write the current set of oligos to file (enables resuming)
        oligo_file = open(file_prefix + 'latest.txt', 'w')
        for i in oligo_list:
            oligo_file.write(i + '\n')
        oligo_file.close()

        # Update the data log file
        data_file = open(file_prefix + 'data.csv', 'a')
        data_csv = csv.writer(data_file, quoting=csv.QUOTE_MINIMAL)
        d = [loop_count, n_attempted, mean_score, min_score, met, timer]
        data_csv.writerow(d)
        data_file.close()


def _n_p_nupack(sequence_list):
    all_concs = []
    for x in sequence_list:
        n_np = Nupack(x, 'dna')
        n_concs = n_np.concentrations(2)['concentration']
        n_concs = n_concs[0:len(x)]
        all_concs.append(n_concs)
    return all_concs


def _monomers_concentration(input):
    sequence_list = input[0]
    conc = input[1]

    # Run nupack's 'concentrations'
    np = Nupack(sequence_list, 'dna')
    nc = np.concentrations(2, conc=conc, mfe=True)
    concs = nc['concentration']
    types = nc['types']

    # Isolate the unbound monomer concentrations
    free_conc = sum([concs[i] for i, x in enumerate(types) if sum(x) == 1])
    free_fraction = free_conc / (2 * conc)
    np.close()  # Delete temp dir

    return free_fraction


def _gen_oligo(prot_seq, min_score, freq_threshold, weighted, conc):
    # Generate oligos until one has a self-self interaction below
    # the threshold.
    threshold_met = False

    while not threshold_met:
        # Generate a randomly or weighted-randomly
        # coding sequence for the specified peptide
        # Ensure it has an mfe of 0
        if weighted:
            choice = weighted_codons.WeightedCodons(prot_seq,
                                                    frequency_table='sc',
                                                    material='pep')
            new_oligo = choice.generate()
        else:
            ft = freq_threshold
            choice = random_codons.RandomCodons(prot_seq,
                                                threshold=ft)
            new_oligo = choice.generate()

        # Check oligo for self-self binding and mfe
        oligo_monomer = _monomers_concentration(([new_oligo, new_oligo], conc))
        oligo_np = Nupack([new_oligo], 'dna')
        oligo_mfe = oligo_np.mfe(0)
        oligo_np.close()
        if oligo_monomer >= min_score and oligo_mfe == 0.0:
            threshold_met = True

    return new_oligo


def _multiprocessing_progress(mp_iterator, total_jobs, interval=1):
    time_start = time.time()
    while (True):
        completed = mp_iterator._index
        if (completed == total_jobs):
            print '\n'
            break
        if completed > 0 & completed % 20 == 0:
            time.sleep(interval)
            remaining = total_jobs - completed
            time_current = time.time()
            rate_current = float(completed) / (time_current - time_start)
            time_remaining = remaining / rate_current
            m = int(time_remaining) / 60  # minutes
            s = time_remaining - 60 * m  # seconds
            print "Remaining: %s:%.2f min(s)\r" % (m, s)
        else:
            time.sleep(.1)


def remove_redundant(seq_list, oligo_params):
    while True:
        counted = collections.Counter(seq_list)
        redundant = []
        new_oligos = []
        for key, value in counted.iteritems():
            if value > 1:
                for i in range(value - 1):
                    redundant.append(key)
        if redundant:
            for x in redundant:
                # Remove redundant entry - the one towards the end of the list
                r_index = len(seq_list) - 1 - seq_list[::-1].index(x)
                seq_list.pop(r_index)
                new_oligos.append(_gen_oligo(*oligo_params))
            seq_list += new_oligos
        else:
            break
    return seq_list


def unflatten(input_list, n):
    list_copy = [x for x in input_list]
    unflattened = []
    for i in range(len(list_copy) / n):
        temp_list = []
        for j in range(n):
            temp_list.append(list_copy.pop(0))
        unflattened.append(temp_list)
    return unflattened


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
    min_score = config.getfloat('Main Settings', 'min unbound monomer score')
    con = config.getfloat('Main Settings', 'species concentration')
    wholemat = config.getboolean('Main Settings', 'whole matrix')
    nm = config.getfloat('Main Settings', 'no increase limit')
    resume = config.get('Main Settings', 'resume')

    # Run OrthoSeq.orthogonal_sequences
    o_seq = OrthoSeq(sequence, T=T, min_score=min_score, n_attempt=nm)
    if resume == 'False':
        o_seq.start(n, m=m, wholemat=wholemat)
    else:
        resfile = open(resume, 'r')
        resume = resfile.readlines()
        resfile.close()
        resume = [x.strip() for x in resume]
        o_seq.start(n, m=m, wholemat=wholemat, resume=resume)
