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

    def __init__(self, seq, dna_conc=50, salt_conc=50,
                 parameters='cloning'):
        '''
        :param seq: Sequence for which to calculate the Tm.
        :type seq: pymbt.sequence.DNA
        :param dna_conc: DNA concentration in nM.
        :type dna_conc: float
        :param salt_conc: Salt concentration in mM.
        :type salt_conc: float
        :param parameters: Nearest-neighbor parameter set. Available options:
                           'breslauer': Breslauer86 parameters
                           'santalucia': SantaLuca98 parameters
                           'cloning': breslauer without corrections
        :type parameters: str

        '''

        self.template = seq

        # store params as attributes for modification?
        self.dna_conc = dna_conc
        self.salt_conc = salt_conc
        self.parameters = parameters

    def run(self):
        '''
        Execute class function.

        '''
        melt = tm(self.template, self.dna_conc, self.salt_conc,
                  self.parameters)
        return melt


def tm(seq, dna_conc=50, salt_conc=50, parameters='cloning'):
    '''
    Returns DNA/DNA tm using nearest neighbor thermodynamics.

    :param seq: Sequence for which to calculate the tm.
    :type seq: pymbt.sequence.DNA
    :param dna_conc: DNA concentration in nM.
    :type dna_conc: float
    :param salt_conc: Salt concentration in mM.
    :type salt_conc: float
    :param parameters: Nearest-neighbor parameter set. Available options:
                       'breslauer': Breslauer86 parameters
                       'santalucia': SantaLuca98 parameters
                       'cloning': breslauer without corrections
    :type parameters: str

    '''

    if parameters == 'cloning':
        params = tm_params.CLONING
    elif parameters == 'breslauer':
        params = tm_params.BRESLAUER
    elif parameters == 'santalucia':
        params = tm_params.SANTALUCIA98
    else:
        raise ValueError('Unsupported parameter set.')

    # Thermodynamic parameters
    pars = {'delta_h': params['delta_h'], 'delta_s': params['delta_s']}
    pars_error = {'delta_h': params['delta_h_err'],
                  'delta_s': params['delta_s_err']}

    # Error corrections - done first for use of reverse_complement parameters
    if parameters == 'breslauer':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'cloning':
        deltas = breslauer_corrections(seq, pars_error)
        deltas[0] += 3.4
        deltas[1] += 12.4
    elif parameters == 'santalucia':
        deltas = santalucia_corrections(seq, pars_error)

    # Sum up the nearest-neighbor enthalpy and entropy
    seq = str(seq).upper()

    def pair_deltas(seq):
        delta0 = 0
        delta1 = 0
        for i in range(len(seq) - 1):
            curchar = seq[i:i + 2]
            delta0 += pars['delta_h'][curchar]
            delta1 += pars['delta_s'][curchar]
        return delta0, delta1

    new_delt = pair_deltas(seq)
    deltas[0] += new_delt[0]
    deltas[1] += new_delt[1]

    # Unit corrections
    salt_conc /= 1e3
    dna_conc /= 1e9
    deltas[0] *= 1e3

    # Universal gas constant (R)
    gas_constant = 1.9872

    if parameters == 'breslauer' or parameters == 'cloning':
        salt_adjusted = 16.6 * log(salt_conc) / log(10.0)
        pre_salt = -deltas[0] / (gas_constant * log(dna_conc / 16.0) -
                                 deltas[1])
        melt = pre_salt + salt_adjusted - 273.15
    else:
        salt_adjusted = 0.368 * (len(seq) - 1) * log(salt_conc) - deltas[1]
        melt = -deltas[0] / (salt_adjusted + gas_constant *
                             log(dna_conc / 4.0)) - 273.15

    return melt


def _pair_count(sequence, pattern):
    '''
    Collect and count sequence pairs in the sequence.

    :param sequence: The string in which to count instances of the pattern.
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


def breslauer_corrections(seq, pars_error):
    deltas_corr = [0, 0]
    contains_gc = 'G' in seq.top or 'C' in seq.top
    only_at = seq.top.count('a') + seq.top.count('t') == len(seq)
    symmetric = seq == seq.reverse_complement()
    terminal_t = seq.top.startswith('t') + seq.top.endswith('t')

    for i, delta in enumerate(['delta_h', 'delta_s']):
        if contains_gc:
            deltas_corr[i] += pars_error[delta]['anyGC']
        if only_at:
            deltas_corr[i] += pars_error[delta]['onlyAT']
        if symmetric:
            deltas_corr[i] += pars_error[delta]['symmetry']
        if terminal_t and delta == 'delta_h':
            deltas_corr[i] += pars_error[delta]['terminalT'] * terminal_t

    return deltas_corr


def santalucia_corrections(seq, pars_error):
    deltas_corr = [0, 0]
    first = seq.top[0]
    last = seq.top[-1]

    startGC = first == 'g' or first == 'c'
    startAT = first == 'a' or first == 't'
    endGC = last == 'g' or last == 'c'
    endAT = last == 'a' or last == 't'
    initGC = startGC + endGC
    initAT = startAT + endAT

    symmetric = seq == seq.reverse_complement()

    for i, delta in enumerate(['delta_h', 'delta_s']):
        deltas_corr[i] += initGC * pars_error[delta]['initGC']
        deltas_corr[i] += initAT * pars_error[delta]['initAT']
        if symmetric:
            deltas_corr[i] += pars_error[delta]['symmetry']

    return deltas_corr
