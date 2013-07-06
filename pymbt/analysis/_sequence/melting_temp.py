'''
Calculate the thermodynamic melting temperatures of nucleotide sequences
using the Finnzymes modified Breslauer 1986 parameters.

'''

from math import log
from pymbt.analysis._sequence import tm_params


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
        melt = tm(str(self.template), self.dna_conc, self.salt_conc,
                  self.method)
        return melt


def tm(sequence, dna_conc=50, salt_conc=50, method='finnzymes'):
    '''
    Returns DNA/DNA tm using nearest neighbor thermodynamics.

    :param sequence: DNA sequence.
    :type sequence: DNA
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
    else:
        msg = "'finnzymes' is the only method that is currently supported"
        raise ValueError(msg)

    # Thermodynamic parameters
    pars = {'delta_h': params['delta_h'], 'delta_s': params['delta_s']}
    pars_error = {'delta_h': params['delta_h_err'],
                  'delta_s': params['delta_s_err']}

    # Error corrections - done first for use of reverse_complement method
    deltas = [0, 0]
    at_count = [x for x in sequence.top if x == 'A' or x == 'T']
    for i, delta in enumerate(['delta_h', 'delta_s']):
        if 'G' in sequence.top or 'C' in sequence.top:
            deltas[i] += pars_error[delta]['initGC']
        if len(at_count) == len(sequence):
            deltas[i] += pars_error[delta]['initAT']
        if sequence.top.startswith('T') and delta == 'delta_h':
            deltas[i] += pars_error[delta]['5termT']
        if sequence == sequence.reverse_complement():
            deltas[i] += pars_error[delta]['symm']

    # Sum up the nearest-neighbor enthalpy and entropy
    sequence = str(sequence).upper()
#    for i, delta in enumerate(['delta_h', 'delta_s']):
#        keys = pars[delta].keys()
#        deltas[i] += sum(_pair_count(sequence, key) * pars[delta][key]
#                         for key in keys)

    # This method is off from the previous one by ~1e-13% from the previous.
    # Not sure why but it doesn't matter. This method is 2X faster.
    def pair_deltas(sequence):
        delta0 = 0
        delta1 = 0
        for i in range(len(sequence) - 1):
            curchar = sequence[i:i + 2]
            delta0 += pars['delta_h'][curchar]
            delta1 += pars['delta_s'][curchar]
        return delta0, delta1

    new_delt = pair_deltas(sequence)
    deltas[0] += new_delt[0]
    deltas[1] += new_delt[1]

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

    # Note: 85% of the tm function's time is spent here - this is the only
    # place to consider optimizing speed

    # Faster than Biopython's overcount
    count = 0
    start = 0
    while True:
        start = sequence.find(pattern, start) + 1
        if start > 0:
            count += 1
        else:
            return count
