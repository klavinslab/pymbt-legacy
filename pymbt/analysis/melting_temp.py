'''Calculate the thermodynamic melting temperatures of nucleotide sequences
using the Finnzymes modified Breslauer 1986 parameters.'''

from math import log
from pymbt.sequence_utils import reverse_complement

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
        params = FINNZYMES_PARAMS
    #elif method == 'santalucia98':
    #    params = SANTALUCIA98_PARAMS
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

FINNZYMES_PARAMS = {
    'delta_h': {
        'AA': 9.1,
        'TT': 9.1,
        'AT': 8.6,
        'TA': 6.0,
        'CA': 5.8,
        'TG': 5.8,
        'GT': 6.5,
        'AC': 6.5,
        'CT': 7.8,
        'AG': 7.8,
        'GA': 5.6,
        'TC': 5.6,
        'CG': 11.9,
        'GC': 11.1,
        'GG': 11.0,
        'CC': 11.0},
    'delta_h_err': {
        'initAT': 0.0,
        'initGC': 0.0,
        'symm': 0.0,
        '5termT': 0.0},
    'delta_s': {
        'AA': 24.0,
        'TT': 24.0,
        'AT': 23.9,
        'TA': 16.9,
        'CA': 12.9,
        'TG': 12.9,
        'GT': 17.3,
        'AC': 17.3,
        'CT': 20.8,
        'AG': 20.8,
        'GA': 13.5,
        'TC': 13.5,
        'CG': 27.8,
        'GC': 26.7,
        'GG': 26.6,
        'CC': 26.6},
    'delta_s_err': {
        'initAT': 0.0,
        'initGC': 0.0,
        'symm': 0.0,
        '5termT': 0.0}}

SANTALUCIA98_PARAMS = {
    'delta_h': {
        'AA': 7.9,
        'TT': 7.9,
        'AT': 7.2,
        'TA': 7.2,
        'CA': 8.5,
        'TG': 8.5,
        'GT': 8.4,
        'AC': 8.4,
        'CT': 7.8,
        'AG': 7.8,
        'GA': 8.2,
        'TC': 8.2,
        'CG': 10.6,
        'GC': 9.8,
        'GG': 8.0,
        'CC': 8.0},
    'delta_h_err': {
        'initAT': 0.0,
        'initGC': 0.0,
        'symm': 0.0,
        '5termT': -0.4},
    'delta_s': {
        'AA': 23.6,
        'TT': 23.6,
        'AT': 18.8,
        'TA': 18.5,
        'CA': 19.3,
        'TG': 19.3,
        'GT': 23.0,
        'AC': 23.0,
        'CT': 16.1,
        'AG': 16.1,
        'GA': 20.3,
        'TC': 20.3,
        'CG': 25.5,
        'GC': 28.4,
        'GG': 15.6,
        'CC': 15.6},
    'delta_s_err': {
        'initAT': 9.0,
        'initGC': 5.9,
        'symm': 1.4,
        '5termT': 0.0}}
