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

import random
import os
import time
import multiprocessing
import socket
from datetime import datetime
from string import maketrans, translate
from itertools import combinations
from math import factorial

#TODO: shouldn't pickle bound methods - use hidden functions (_function)
import copy_reg
import types  # necessary to pickle bound methods

from pymbt.nupack import Nupack
from pymbt.dna_manipulation import reverse_complement as r_c

# To do:
# 1. Implement signal handler so it exits cleanly on ctrl-c
# 2. Group NxM and MxM calculations before parallelizing
# 3. Many parts of the while loop are repetitive
#        -roll into for loop or something.
# 4. self-self is calculated twice right now - remove that.


class OrthoSeq:
    def __init__(self,
                 prot_seq,
                 T=50,
                 min_free=0.95,
                 conc=5e-7,
                 nullinc_max=1e5):
        # Make sure sequence only includes amino acids
        protein_alphabet = "ACDEFGHIKLMNPQRSTVWY"
        for i, char in enumerate(prot_seq):
            if char not in protein_alphabet:
                raise ValueError("Non-amino acid character at index " + str(i))

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
                             report_nupack=False):

        # The default combination number should be done dynamically.
        # Solved combinations equation for m.
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
            n_oligos = [x for x in o_iterator]
            n_pool.close()
            n_pool.join()

            # Now calculate larger pool of oligos.
            # It is probably safe to only calculate:
            #   1. Each oligo in a 'forward' version against all
            #      other reversed oligos in a size-N pool
            #   2. Repeat this for each oligo
            #   3. Discard the m worst oligos, repeat

            oligo_list += new_oligos

            # Generate all the combinations of oligos in this set
            combos = [list(x) for x in combinations(oligo_list, n)]

            # For each combination of length n,
            # make n lists of oligos where each one
            # is paired with the reverse of all the others
            full_combos = [[] for i, x in enumerate(combos)]
            for i, v in enumerate(combos):
                for j in v:
                    newlist = [j]
                    newlist.append([rc(x) for x in v if x is not j])
                    full_combos[i].append(newlist)

            # The list generated is of length n+1 containing lists of length n
            # In order to feed to nupack, have to flatten the list
            # this may only be true for m = 1
            full_combos = [v for x in full_combos for v in x]

            # Run that list through nupack
            o_pool = multiprocessing.Pool()
            o_pool_iter = o_pool.imap(self.n_p_nupack, full_combos)
            # Report progress
            if report_nupack or n + m > 6:
                tot = len(full_combos)
                self._multiprocessing_progress(o_pool_iter, tot, interval=5)
            mon_concs = [x for x in o_pool_iter]
            o_pool.close()
            o_pool.join()

            # Group results into n+1 lists of length n lists
            np_result = []
            while len(mon_concs) is not 0:
                np_result.append([mon_concs.pop(0) for i in range(n)])
            np_result_flat = [[w for v in x for w in v] for x in np_result]

            # Keep the best-scoring combo
            mc_min = [min(x) for x in np_result_flat]
            best_combo = mc_min.index(max(mc_min))
            oligo_list = combos[best_combo]
            mon_concs_left = np_result_flat[best_combo]

            # Remove the one with the worst score
            #mc = [sum(x) for x in mon_concs]
            #for i in range(m):
            #    oligo_list.pop(mc.index(min(mc)))
            #    mon_concs.pop(mc.index(min(mc)))
            #    mc.pop(mc.index(min(mc)))

            ###########
            # Scoring
            ###########
            mean_free = (sum(mon_concs_left) / len(mon_concs_left)) / self.conc
            min_free = min(mon_concs_left) / self.conc

            met = 0
            for x in np_result[best_combo]:
                if min(x) / self.conc >= self.min_free:
                    met += 1
            if min_free == last_min_free:
                nullinc += 1
            else:
                nullinc = 0

            #######################################
            # Logging - provides data trail and allows resuming/passing along
            # an extended attempt
            #######################################

            # Initial setup steps
            if loop_count == 1:
                cur_date = datetime.now().strftime('%Y%m%d--%H%M%S')
                cur_dir = os.getcwd()
                host = socket.gethostname()
                cur_id = cur_date + '-' + host + '-' + self.prot_seq
                fileprefix = cur_dir + '/' + cur_id + '-oligo_calc_'

                # Write out setup info file
                infofile = open(fileprefix + 'info.txt', 'w')
                infofile.write('sequence oligo_n oligo_m \n')
                infofile.write(self.prot_seq + ' %i %i \n' % (n, m))
                infofile.close()

                # Set up data log file
                datafile = open(fileprefix + 'data.txt', 'w')
                header = 'loop nullinc mean_free min_free met time\n'
                datafile.write(header)

            # Write the current set of oligos to file (enables resuming)
            oligofile = open(fileprefix + 'latest.txt', 'w')
            for i in oligo_list:
                oligofile.write(i + '\n')
            oligofile.close()

            # Update the data log file
            datafile = open(fileprefix + 'data.txt', 'a')
            time_diff = time.time() - loop_start
            time_elapsed += time_diff
            d = [loop_count, nullinc, mean_free, min_free, met, time_elapsed]
            datafile.write(' '.join(d))
            datafile.write('\n')
            datafile.close()

            # Stop looping if the threshold is met
            if met == len(oligo_list):
                break
            if nullinc == self.nullinc_max:
                break

            last_min_free = min_free
            """

            #######################################
            # Script Stop decision
            #######################################
            # Stop when all oligos interact less than min threshold
            if worst_free > self.min_free:
                break

            # If haven't increased score in 'noinc' turns, break

            if stop == loop_count:
                break
            """
        return oligo_list

    def n_p_nupack(self, sequence_list):
        n_np = Nupack(sequence_list, 'dna')
        n_concs = n_np.concentrations(2)['concentration'][0:len(sequence_list)]
        return n_concs

    def gen_oligo(self, min_free, min_freq=0.3, weighted=False):
        # Dict of dicts to hold codon usage information.
        # Source of codon frequencies:
        #   http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932
        c_dict = {
            'A': {'GCG': 0.109972396541529,
                  'GCA': 0.288596474496094,
                  'GCT': 0.377014739102356,
                  'GCC': 0.224416389860021},
            'R': {'AGG': 0.208564104515562,
                  'AGA': 0.481137590939125,
                  'CGG': 0.0392677130215486,
                  'CGA': 0.0676728924436203,
                  'CGT': 0.144572019635586,
                  'CGC': 0.0587856794445578},
            'N': {'AAT': 0.589705127199784,
                  'AAC': 0.410294872800217},
            'D': {'GAT': 0.65037901553924,
                  'GAC': 0.34962098446076},
            'C': {'TGT': 0.629812614586062,
                  'TGC': 0.370187385413938},
            '*': {'TGA': 0.303094329334787,
                  'TAG': 0.225736095965104,
                  'TAA': 0.471169574700109},
            'Q': {'CAG': 0.307418833439535,
                  'CAA': 0.692581166560465},
            'E': {'GAG': 0.296739610207218,
                  'GAA': 0.703260389792782},
            'G': {'GGG': 0.119057918187951,
                  'GGA': 0.215422869017838,
                  'GGT': 0.472217600813099,
                  'GGC': 0.193301611981112},
            'H': {'CAT': 0.636710255236351,
                  'CAC': 0.363289744763649},
            'I': {'ATA': 0.273331091899568,
                  'ATT': 0.462925823433014,
                  'ATC': 0.263743084667417},
            'L': {'TTG': 0.286319859527146,
                  'TTA': 0.275534472444779,
                  'CTG': 0.110440170850593,
                  'CTA': 0.141277445174148,
                  'CTT': 0.129115062940288,
                  'CTC': 0.0573129890630467},
            'K': {'AAG': 0.423936637198697,
                  'AAA': 0.576063362801303},
            'M': {'ATG': 1},
            'F': {'TTT': 0.586126603840976,
                  'TTC': 0.413873396159024},
            'P': {'CCG': 0.120626895854398,
                  'CCA': 0.417143753704543,
                  'CCT': 0.307740315888567,
                  'CCC': 0.154489034552491},
            'S': {'AGT': 0.159245398699046,
                  'AGC': 0.109749229743856,
                  'TCG': 0.0963590866114069,
                  'TCA': 0.210157220085731,
                  'TCT': 0.264456618519558,
                  'TCC': 0.160032446340401},
            'T': {'ACG': 0.135583991997041,
                  'ACA': 0.302413913478422,
                  'ACT': 0.345237040780705,
                  'ACC': 0.216765053743832},
            'W': {'TGG': 1},
            'Y': {'TAT': 0.559573963633711,
                  'TAC': 0.440426036366289},
            'V': {'GTG': 0.190897642582249,
                  'GTA': 0.208783185960798,
                  'GTT': 0.391481704636128,
                  'GTC': 0.208837466820824}}

        # Remove codons that don't meet the threshold
        cd = {}
        for i in c_dict:
            max_i = max(c_dict[i].values())
            cd[i] = [x for x, v in c_dict[i].items() if v / max_i >= min_freq]

        # Ensure that trimmed dictionary includes the amino acids needed
        for i, char in enumerate(self.prot_seq):
            if len(cd[char]) == 0:
                raise ValueError("No codons meet the frequency threshold " +
                "for residue " + str(char))

        # Generate oligos until one has a self-self interaction below
        # the threshold.
        # There is currently no timeout set - could go on forever if
        # min_free is too stringent
        oligo_meets_free_conc = False

        while not oligo_meets_free_conc:
            # Generate a random oligo that
            # has an mfe score of 0 (no self-folding)
            if weighted:
                def weightedchoice(aminoacid):
                    # cd contains allowed codons
                    # c_dict contains usage frequencies
                    simp = [[x, c_dict[aminoacid][x]] for x in cd[aminoacid]]
                    total = sum([x[1] for x in simp])
                    rn = random.uniform(0, total)
                    cumsum = 0
                    for x, v in simp:
                        if cumsum + v > rn:
                            return x
                        else:
                            cumsum += v
                n_list = [weightedchoice(x) for x in self.prot_seq]
                new_oligo = ''.join(n_list)
            else:
                n_list = [random.sample(cd[x], 1)[0] for x in self.prot_seq]
                new_oligo = ''.join(n_list)

            # Check oligo self-self binding
            n_oligo = [new_oligo for x in range(2)]
            oligo_monomer = self.monomers_concentration(n_oligo)
            oligo_np = Nupack([new_oligo], 'dna')
            oligo_mfe = oligo_np.mfe(0)
            oligo_np.close()
            if oligo_monomer >= min_free and oligo_mfe == 0.0:
                oligo_meets_free_conc = True

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
        sl = seq_list  # shorter name
        # List generated is actually 1-dimensional, but results are parsed into
        # a list-of-lists matrix format later
        seq_pairs = []
        if reverse:
            for i, x in enumerate(sl):
                for j in range(i, len(sl)):
                    if i is not j:
                        seq_pairs.append([sl[i], rc(sl[j])])
        # seq_pairs = [[sl[i], rc(sl[j])] for i, x
        # in enumerate(sl) for j in range(i, len(sl)) if i is not j]
        else:
            for i, x in enumerate(sl):
                for j in range(i, len(sl)):
                    seq_pairs.append([sl[i], sl[j]])
        # seq_pairs = [[sl[i], sl[j] ] for i, x
        # in enumerate(sl) for j in range(i, len(sl))]

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
                s = time_remaining - 60 * time_mins  # seconds
                print("Remaining: %(m).f:%(s) min(s)\r" % {'m': m, 's': s})
                time.sleep(interval)
            else:
                time.sleep(.1)

    def _score_by_mean_free(self, mat):
        mat_flat = [x for v in mat for x in v if type(x) is not list]
        mat_mean = sum(mat_flat) / len(mat_flat)
        return mat_mean

    def _score_by_min_free(self, mat):
        mat_flat = [x for v in mat for x in v]
        mat_min = min(mat_flat)
        return mat_min


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
