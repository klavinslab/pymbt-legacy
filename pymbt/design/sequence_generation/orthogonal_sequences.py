'''
Generates a set of orthogonal DNA sequences for a given protein sequence
Currently limited by the size of the protein sequence - ~10 residues max
To increase size, need to optimize search method and take steps to reduce
memory usage (write to files and iterate/queue).

'''

import csv
import collections
from datetime import datetime
from itertools import combinations
import os
import socket
import time
from pymbt.analysis import Nupack, nupack_multiprocessing
from pymbt.design import RandomCodons
from pymbt.design import WeightedCodons
from pymbt.sequence.utils import reverse_complement as r_c
from pymbt.sequence.utils import check_alphabet

# TODO: allow OrthoSeq to take config file and previous run as keyword inputs.
# then the config file can be placed in a separate dir

# TODO: make this compatible with DNA objects


class OrthoSeq(object):
    '''Class to calculate, store, and write orthogonality-optimized peptide
    sequences.'''

    def __init__(self, prot_seq, oligo_n, candidates=0, temp=50,
                 min_score=0.95, n_attempt=1e5):
        '''
        :param prot_seq: Input protein sequence.
        :type prot_seq: str
        :param oligo_n: The final number of orthogonal sequences to output.
        :type oligo_n: int
        :param candidates: The number of new sequences to try out on each
                           iteration.
        :type candidates: int
        :param temp: Temperature at which to evaluate for orthogonality.
        :type temp: float
        :param min_score: Score setpoint - stop optimization once this score
                          has been reched.
        :type min_score: float
        :param n_attempt: Attempt timeout - sets maximum number of times to
                          attempt optimization without seeing an improvement in
                          score

        '''

        # Make sure sequence only includes amino acids
        check_alphabet(prot_seq, material='pep')

        # Attributes
        self.prot_seq = prot_seq
        self.oligo_n = oligo_n
        self.candidates = candidates
        self.temp = temp
        # Primary scoring metric - minimum average free monomer concentrations
        # - 0 would be fully-bound, 1 would be nothing bound
        self.min_score = min_score
        # Secondary scoring metric - stops the simulation after n
        # attempts that fail to increase the primary score
        self.n_attempt = n_attempt

    def start(self, resume=False, report=True, weighted=True,
              freq_threshold=0.5, oligo_conc=5e-7):
        '''
        Start the optimization.

        :param resume: A path to a dir from a previous run, and will
                       cause start() to pick up where it left off.
        :type resume: str
        :param report: report multiprocessing progress.
        :type report: bool
        :param weighted: Determines whether codons generated should be
                         weighted (codon-optimized).
        :type weighted: bool
        :param freq_threshold: relative frequency threshold below which codons
                               should not be used. Enables avoidance of very
                               rare codons (codon optimization).
        :type freq_threshold: float

        '''

        # Arguments for _gen_oligo sequence generation
        oligo_args = (self.prot_seq,
                      self.min_score,
                      freq_threshold,
                      weighted,
                      oligo_conc)
        # Generate n oligos that meet score threshold
        # If there's a list of oligos to resume improving, start using it
        if resume:
            if len(resume) != self.oligo_n:
                raise Exception('Number of sequences to resume doesn\'t match\
                                 setting of \'n\'.')
            for oligo in resume:
                if type(oligo) != str:
                    raise Exception('\'resume\' input isn\'t strings!')
            oligo_list = resume
        # If not, generate initial set of oligos
        else:
            oligo_list = [_gen_oligo(*oligo_args) for x in range(self.oligo_n)]
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
            new_oligos = [_gen_oligo(*oligo_args) for x in
                          range(self.candidates)]
            oligo_list += new_oligos
            # If two or more oligos are the same, replace with unique ones
            oligo_list = remove_redundant(oligo_list, oligo_args)
            # Generate all combinations - (m+n) choose n
            combos = [list(x) for x in combinations(oligo_list, self.oligo_n)]

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
                                                  {'max_complexes': 2,
                                                   'report': report})
            # just the concentrations list from each entry
            np_concs = [x['concentration'] for x in np_concs_raw]
            # just the types list from each entry
            np_types = [x['types'] for x in np_concs_raw]
            # the monomer concentrations from each entry
            mon_concs = []
            for np_conc, np_type in zip(np_concs, np_types):
                # if sum of the list is 1, that means it's a monomer conc
                # (e.g. [1, 0, 0, 0])
                temp_list = []
                for j, vals in enumerate(np_type):
                    if sum(vals) == 1:
                        temp_list.append(np_conc[j])
                mon_concs.append(temp_list)

            # Group by combo
            grouped_mon_concs = unflatten(mon_concs, self.candidates)

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
            score_mean = (sum(best_concs) / len(best_concs)) / oligo_conc
            score_min = min(best_concs) / oligo_conc
            met = 0
            for conc in best_concs:
                if (conc / oligo_conc) >= self.min_score:
                    met += 1

            # Did the minimum score improve? If not, increment counter
            if score_min == last_min_score:
                n_attempted += 1
            else:
                n_attempted = 0
                oligo_list = best_combo

            # Record how much time has passed
            timer = time.time() - time_start

            # Log
            self._log(loop_count,
                      n_attempted,
                      score_mean,
                      score_mean,
                      met,
                      timer,
                      oligo_list)

            # Script stop decision:
            # Stop looping if either:
            #   the threshold is met
            #   the minimum score has stopped improving by n_attempt
            if met == len(oligo_list):
                break
            if n_attempted == self.n_attempt:
                break

            last_min_score = score_min

        return oligo_list

    def _log(self, loop_count, n_attempted, score_mean, score_min, met, timer,
             oligo_list):
        '''
        Provides data trail and allows resuming/passing along an extended
        attempt.

        :param loop_count: Current loop number.
        :type loop_count: int
        :param n_attempted: Number of times optimization has been consecutively
                            attempted without improving the score.
        :type n_attempted: int
        :param score_mean: Mean of the unboundedness score.
        :type score_mean: float
        :param score_min: Min of the unboundedness score.
        :type score_min: float
        :param met: Number of oligos that have met the score threshold.
        :type met: int
        :param timer: Length of time for the calculation so far in seconds.
        :type timer: float
        :param oligo_list: List of the current best oligos.
        :type oligo_list: list

        '''

        # Initial file setup for logging
        if loop_count == 1:
            # Use date, time, host, and peptide sequence to uniquify file
            current_date = datetime.now().strftime('%Y%m%d--%H%M%S')
            host = socket.gethostname()
            cur_id = current_date + '-' + host + '-' + self.prot_seq
            cur_dir = os.getcwd()
            file_prefix = '{0}/{1}-oligo_calc_'.format(cur_dir, cur_id)

            # Write out setup info file
            info_file = open(file_prefix + 'info.txt', 'w')
            info_csv = csv.writer(info_file, quoting=csv.QUOTE_MINIMAL)
            info_csv.writerow(['sequence', 'oligo_n', 'candidates'])
            info_csv.writerow([self.prot_seq, self.oligo_n, self.candidates])
            info_file.close()

            # Set up data log file
            data_file = open(file_prefix + 'data.csv', 'w')
            data_csv = csv.writer(data_file, quoting=csv.QUOTE_MINIMAL)
            cols = ['loop',
                    'n_attempted',
                    'score_mean',
                    'score_min',
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
        data = [loop_count, n_attempted, score_mean, score_min, met, timer]
        data_csv.writerow(data)
        data_file.close()


def _n_p_nupack(sequence_list):
    '''
    Runs nupack modules 'concentrations' on a list of sequences. Defined at top
    level to enable pickling for multiprocessing.

    :param sequence_list: Sequences to evaluate.
    :type sequence_list: list

    '''

    all_concs = []
    for sequence in sequence_list:
        n_np = Nupack(sequence, 'dna')
        n_concs = n_np.concentrations(2)['concentration']
        n_concs = n_concs[0:len(sequence)]
        all_concs.append(n_concs)
    return all_concs


def _monomers_concentration(sequences, conc):
    '''
    Calculates the unbound (monomer) concentrations for a given concentration
    of oligos.

    :param input: Tuple of a list of sequences and their concentrations,
                  respectively.
    :type input: tuple

    '''

    # Run nupack's 'concentrations'
    nupack = Nupack(sequences, 'dna')
    concentrations = nupack.concentrations(2, conc=conc, mfe=True)
    concs = concentrations['concentration']
    types = concentrations['types']

    # Isolate the unbound monomer concentrations
    free_conc = sum([c for c, x in zip(concs, types) if sum(x) == 1])
    free_fraction = free_conc / (2 * conc)
    nupack.close()  # Delete temp dir

    return free_fraction


def _gen_oligo(prot_seq, min_score, freq_threshold, weighted, conc):
    '''
    Generate oligos until one has a self-self interaction below the threshold.

    :param prot_seq: Protein sequence.
    :type prot_seq: str
    :param min_score: Score setpoint - stop optimization once this score
                      has been reched.
    :type min_score: float
    :param freq_threshold: relative frequency threshold below which codons
                           should not be used. Enables avoidance of very rare
                           codons (codon optimization).
    :type freq_threshold: float
    :param weighted: Determines whether codons generated should be
                     weighted (codon-optimized).
    :type weighted: bool

    '''

    threshold_met = False

    while not threshold_met:
        # Generate a randomly or weighted-randomly
        # coding sequence for the specified peptide
        # Ensure it has an mfe of 0
        if weighted:
            choice = WeightedCodons(prot_seq,
                                    frequency_table='sc',
                                    material='pep')
            new_oligo = choice.generate()
        else:
            choice = RandomCodons(prot_seq,
                                  threshold=freq_threshold)
            new_oligo = choice.generate()

        # Check oligo for self-self binding and mfe
        oligo_monomer = _monomers_concentration([new_oligo, new_oligo], conc)
        oligo_np = Nupack([new_oligo], 'dna')
        oligo_mfe = oligo_np.mfe(0)
        oligo_np.close()
        if oligo_monomer >= min_score and oligo_mfe == 0.0:
            threshold_met = True

    return new_oligo


def remove_redundant(seq_list, oligo_params):
    '''
    Removes and regenerates sequences that appear more than once in a list.

    :param seq_list: Sequences.
    :type seq_list: str
    :param oligo_params: parameters to be passed to _gen_oligo.
    :type oligo_params: tuple

    '''

    while True:
        counted = collections.Counter(seq_list)
        redundant = []
        new_oligos = []
        for key, value in counted.iteritems():
            if value > 1:
                for i in range(value - 1):
                    redundant.append(key)
        if redundant:
            for sequence in redundant:
                # Remove redundant entry - the one towards the end of the list
                r_index = len(seq_list) - 1 - seq_list[::-1].index(sequence)
                seq_list.pop(r_index)
                new_oligos.append(_gen_oligo(*oligo_params))
            seq_list += new_oligos
        else:
            break
    return seq_list


def unflatten(input_list, size):
    '''
    Given a flat list, makes a new list of 'size'-length lists out of it, in
    order.

    :param input_list: List to flatten.
    :type input_list: list
    :param size: Group size.
    :type size: int

    '''

    list_copy = [x for x in input_list]
    unflattened = []
    for i in range(len(list_copy) / size):
        temp_list = []
        for j in range(size):
            temp_list.append(list_copy.pop(0))
        unflattened.append(temp_list)
    return unflattened


if __name__ == "__main__":
    import ConfigParser

    # Read config
    config = ConfigParser.RawConfigParser()
    config.read(os.path.abspath(os.path.dirname(__file__) +
                '/orthogonal_sequences.cfg'))

    # Parse config into parameters
    S = config.get('Main Settings', 'protein sequence')
    N = config.getint('Main Settings', 'n sequences')
    M = config.getint('Main Settings', 'm sequences')
    T = config.getfloat('Main Settings', 'temperature')
    score_thresh = config.getfloat('Main Settings',
                                   'min unbound monomer score')
    CONC = config.getfloat('Main Settings', 'species concentration')
    NM = config.getfloat('Main Settings', 'no increase limit')
    RESUME = config.get('Main Settings', 'resume')

    # Run OrthoSeq.orthogonal_sequences
    NEW_RUN = OrthoSeq(S, N, candidates=M, T=T, min_score=score_thresh,
                       n_attempt=NM)
    if not RESUME:
        NEW_RUN.start()
    else:
        with open(RESUME, 'r') as resfile:
            PREV = [oligos.strip() for oligos in resfile.readlines()]
        resfile.close()

        NEW_RUN.start(resume=PREV)
