'''Calculate the thermodynamic melting temperatures of nucleotide sequences
using the Finnzymes modified Breslauer 1986 parameters.'''

from math import log
from pymbt.sequence.utils import reverse_complement
from pymbt.analysis.sequence import tm_params

# TODO: Add Breslauer, SantaLucia98, and Sugimoto methods
# See doi: 10.1093/bioinformatics/bti066 for good comparison
#   "Comparison of different melting temperature calculation
#    methods for short DNA sequences"


class Tm(object):
    '''
    Nearest-neighbor melting temperature (Tm) calculator.
    '''
    def __init__(self, dna_object, dna_conc=50, salt_conc=50,
                 method='finnzymes'):
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

        self.template = dna_object

        # store params as attributes for modification?
        self.dna_conc = dna_conc
        self.salt_conc = salt_conc
        self.method = method

    def run(self):
        '''
        Execute class function.

        '''
        melt = _calc_tm(str(self.template), self.dna_conc, self.salt_conc,
                        self.method)
        return melt


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

    sequence = sequence.upper()

    if method == 'finnzymes':
        params = tm_params.FINNZYMES_PARAMS
    else:
        msg = '\'finnzymes\' is the only method that is currently supported'
        raise ValueError(msg)

    # Thermodynamic parameters
    pars = {'delta_h': params['delta_h'], 'delta_s': params['delta_s']}
    pars_error = {'delta_h': params['delta_h_err'],
                  'delta_s': params['delta_s_err']}

    # Sum up the nearest-neighbor enthalpy and entropy
    deltas = []
    for delta in ['delta_h', 'delta_s']:
        keys = pars[delta].keys()
        deltas.append(sum(_pair_count(sequence, key) * pars[delta][key]
                      for key in keys))

    # Error corrections
    at_count = [x for x in sequence if x == 'A' or x == 'T']
    for i, delta in enumerate(['delta_h', 'delta_s']):
        if 'G' in sequence or 'C' in sequence:
            deltas[i] += pars_error[delta]['initGC']
        if len(at_count) == len(sequence):
            deltas[i] += pars_error[delta]['initAT']
        if sequence.startswith('T') and delta == 'delta_h':
            deltas[i] += pars_error[delta]['5termT']
        if sequence == reverse_complement(sequence, 'dna'):
            deltas[i] += pars_error[delta]['symm']

    # Unit corrections
    salt_conc /= 1e3
    dna_conc /= 1e9
    deltas[0] *= 1e3

    salt_conc_adjusted = 16.6 * log(salt_conc) / log(10.0)

    # Universal gas constant (R)
    gas_constant = 1.9872

    # These corrections are unaccounted for but are required
    # for the 'finnzymes' method to be accurate.
    if method == 'finnzymes':
        deltas[0] += 3400
        deltas[1] += 12.4
        pre_salt = -deltas[0] / (gas_constant * log(dna_conc / 16) - deltas[1])
        melt = pre_salt + salt_conc_adjusted - 273.15

    return melt


def _pair_count(sequence, pattern):
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
