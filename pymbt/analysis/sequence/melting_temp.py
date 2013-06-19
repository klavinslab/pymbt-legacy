'''Calculate the thermodynamic melting temperatures of nucleotide sequences
using the Finnzymes modified Breslauer 1986 parameters.'''

from math import log
from pymbt.sequence.utils import reverse_complement, check_instance
from pymbt.analysis.sequence import tm_params

# TODO: Add Breslauer, SantaLucia98, and Sugimoto methods
# See doi: 10.1093/bioinformatics/bti066 for good comparison
#   "Comparison of different melting temperature calculation
#    methods for short DNA sequences"


class Tm(object):
    '''
    :param dna_object: DNA sequence.
    :type dna_object: DNA object
    :param dna_conc: DNA concentration in nM.
    :type dna_conc: float
    :param salt_conc: Salt concentration in mM.
    :type salt_conc: float
    :param method: computation method to use. Only available method is
                   'finnzymes'.
    :type method: str
    '''
    def __init__(self, dna_object, dna_conc=50, salt_conc=50,
                 method='finnzymes'):
        # TODO: type checking / dsDNA checking
        self.template = dna_object
        check_instance(self.template)

        # store params as attributes for modification?
        self.dna_conc = dna_conc
        self.salt_conc = salt_conc
        self.method = method

    def run(self):
        tm = _calc_tm(str(self.template), self.dna_conc, self.salt_conc,
                      self.method)
        return tm


def _calc_tm(sequence, dna_conc=50, salt_conc=50, method='finnzymes'):
    '''
    Returns DNA/DNA tm using nearest neighbor thermodynamics.

    :param sequence: DNA sequence.
    :type sequence: str
    :param dna_conc: DNA concentration in nM.
    :type dna_conc: float
    :param salt_conc: Salt concentration in mM.
    :type salt_conc: float
    :param method: computation method to use. Only available method is
                   'finnzymes'.
    :type method: str

    '''

    if method == 'finnzymes':
        params = tm_params.FINNZYMES_PARAMS
    #elif method == 'santalucia98':
    #    params = tm_params.SANTALUCIA98_PARAMS
    else:
        msg = '\'finnzymes\' is the only method that is currently supported'
        raise ValueError(msg)
    delta_h_par = params['delta_h']
    delta_s_par = params['delta_s']
    delta_h_par_err = params['delta_h_err']
    delta_s_par_err = params['delta_s_err']

    def collect_pairs(sequence, pattern):
        '''
        Collect and count sequence pairs in the sequence.

        :param sequence: Any string (in this case, DNA or RNA).
        :type sequence: str
        :param pattern: Pair of strings to search for (in this case,
                        two bases).
        :type pattern: str

        '''

        # Faster than Biopython's overcount
        count = 0
        start = 0
        while True:
            start = sequence.find(pattern, start) + 1
            if start > 0:
                count += 1
            else:
                return count

    # shorter name
    seq = sequence
    seq = seq.upper()
    # Sum up the nearest-neighbor enthalpy and entropy counts*par
    h_keys = delta_h_par.keys()
    s_keys = delta_s_par.keys()
    delta_h = sum([collect_pairs(seq, x) * delta_h_par[x] for x in h_keys])
    delta_s = sum([collect_pairs(seq, x) * delta_s_par[x] for x in s_keys])

    # Error corrections
    # Initiation for seqs with a G-C pair vs only A-T
    if 'G' in seq or 'C' in seq:
        delta_h += delta_h_par_err['initGC']
        delta_s += delta_s_par_err['initGC']
    at_count = [x for x in seq if x == 'A' or x == 'T']
    if len(at_count) == len(seq):
        delta_h += delta_h_par_err['initAT']
        delta_s += delta_s_par_err['initAT']
    # 5' terminal T-A
    if seq.startswith('T'):
        delta_h += delta_h_par_err['5termT']
        delta_s += delta_s_par_err['5termT']
    # Correction for self-complementary sequences
    # The meaning of 'self-complementary' is not well-defined...
    if seq == reverse_complement(seq):
        delta_h += delta_h_par_err['symm']
        delta_s += delta_s_par_err['symm']

    # Unit corrections
    salt_conc = salt_conc / 1e3
    salt_conc_adjusted = 16.6 * log(salt_conc) / log(10.0)
    dna_conc = dna_conc / 1e9
    delta_h = delta_h * 1e3

    # Universal gas constant (R)
    gas_constant = 1.9872

    # These corrections are unaccounted for but are required
    # for the 'finnzymes' method to be accurate.
    if method == 'finnzymes':
        delta_h += 3400
        delta_s += 12.4
        pre_salt = -delta_h / (gas_constant * log(dna_conc / 16) - delta_s)
        melting_temp = pre_salt + salt_conc_adjusted - 273.15
    if method == 'santalucia98':
        k = dna_conc / 4.0
        delta_s = delta_s - 0.368 * (len(seq) - 1) * log(salt_conc)
        melting_temp = -delta_h / ((gas_constant * log(k)) - delta_s) - 273.15

    return melting_temp
