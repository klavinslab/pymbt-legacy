'''
Calculate the thermodynamic melting temperatures of nucleotide sequences.

'''

# TODO: Owczarzy et al 2004 has better salt correction
# TODO: Remove sugimoto? It's missing important details (like salt correction)
# TODO: Write adjusted santalucia method. Did a fit of 20000 sequences
# comparing santalucia98 to finnzymes' modified breslauer method. As expected
# they correlate heavily with offset of -6.014 and slope of 1.335 -
# i.e. to get approximately the same result, do santalucia98, multiply by 1.335
# , and subtract 6.
# Double check those stats
# TODO: Make new hybrid method - combine santalucia unified with owczarzy
# corrections, compare to finnzymes.

from math import log, log10
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
                           'sugimoto': Sugimoto96 parameters
                           'santalucia96': SantaLucia96 parameters
                           'santalucia98': SantaLucia98 parameters
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
                       'sugimoto': Sugimoto96 parameters
                       'santalucia96': SantaLucia96 parameters
                       'santalucia98': SantaLucia98 parameters
                       'cloning': breslauer without corrections
    :type parameters: str

    '''

    if parameters == 'breslauer':
        params = tm_params.BRESLAUER
    elif parameters == 'sugimoto':
        params = tm_params.SUGIMOTO
    elif parameters == 'santalucia96':
        params = tm_params.SANTALUCIA96
    elif parameters == 'santalucia98':
        params = tm_params.SANTALUCIA98
    elif parameters == 'cloning':
        params = tm_params.CLONING
    else:
        raise ValueError('Unsupported parameter set.')

    # Thermodynamic parameters
    pars = {'delta_h': params['delta_h'], 'delta_s': params['delta_s']}
    pars_error = {'delta_h': params['delta_h_err'],
                  'delta_s': params['delta_s_err']}

    # Error corrections - done first for use of reverse_complement parameters
    if parameters == 'breslauer':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'sugimoto':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'santalucia96':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'santalucia98':
        deltas = santalucia98_corrections(seq, pars_error)
    elif parameters == 'cloning':
        deltas = breslauer_corrections(seq, pars_error)
        deltas[0] += 3.4
        deltas[1] += 12.4

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
    R = 1.9872

    # Supposedly this is what dnamate does, but the output doesn't match theirs
#    melt = (-deltas[0] / (-deltas[1] + R * log(dna_conc / 4.0))) +
#                          16.6 * log(salt_conc) - 273.15
#    return melt
    # Overall equation is supposedly:
    # sum{dH}/(sum{dS} + R ln(dna_conc/b)) - 273.15
    # with salt corrections for the whole term (or for santalucia98,
    # salt corrections added to the dS term.
    # So far, implementing this as described does not give results that match
    # any calculator but Biopython's

    if parameters == 'breslauer' or parameters == 'cloning':
        numerator = -deltas[0]
        # Modified dna_conc denominator
        denominator = (-deltas[1]) + R * log(dna_conc / 16.0)
        # Modified Schildkraut-Lifson equation adjustment
        salt_adjustment = 16.6 * log(salt_conc) / log(10.0)
        melt = numerator / denominator + salt_adjustment - 273.15
    elif parameters == 'santalucia98':
        # This one is definitely correct
        # TODO: dna_conc should be divided by 2.0 when dna_conc >> template
        # (like PCR)
        numerator = -deltas[0]
        # SantaLucia 98 salt correction
        salt_adjustment = 0.368 * (len(seq) - 1) * log(salt_conc)
        denominator = -deltas[1] + salt_adjustment + R * log(dna_conc / 4.0)
        melt = -deltas[0] / denominator - 273.15
    elif parameters == 'santalucia96':
        # TODO: find a way to test whether the code below matches another
        # algorithm. It appears to be correct, but need to test it.
        numerator = -deltas[0]
        denominator = -deltas[1] + R * log(dna_conc / 4.0)
        # SantaLucia 96 salt correction
        salt_adjustment = 12.5 * log10(salt_conc)
        melt = numerator / denominator + salt_adjustment - 273.15
    elif parameters == 'sugimoto':
        # TODO: the stuff below is untested and probably wrong
        numerator = -deltas[0]
        denominator = -deltas[1] + R * log(dna_conc / 4.0)
        # Sugimoto parameters were fit holding salt concentration constant
        # Salt correction can be chosen / ignored? Remove sugimoto set since
        # it's so similar to santalucia98?
        salt_correction = 16.6 * log10(salt_conc)
        melt = numerator / denominator + salt_correction - 273.15

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
    terminal_t = seq.top[0] == 't' + seq.top[-1] == 't'

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


def santalucia98_corrections(seq, pars_error):
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
