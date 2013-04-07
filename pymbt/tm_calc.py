# TODO: Add Breslauer, SantaLucia98, and Sugimoto methods
# TODO: Note: santalucia98 is broken
# TODO: Redo santalucia NN parameters. Find a standard against which to compare
# TODO: review unit corrections, see which apply for finnzymes vs. others
# See doi: 10.1093/bioinformatics/bti066 for good comparison
#   "Comparison of different melting temperature calculation
#    methods for short DNA sequences"

'''Calculate the thermodynamic melting temperatures of nucleotide sequences
using the Finnzymes modified Breslauer 1986 parameters.'''

from math import log
from pymbt.sequence_manipulation import reverse_complement


def calc_tm(sequence, dnac=50, saltc=50, method='finnzymes'):
    '''
    Returns DNA/DNA tm using nearest neighbor thermodynamics.

    :param sequence: DNA sequence.
    :type sequence: str.
    :param dnac: DNA concentration in nM.
    :type dnac: float.
    :param saltc: Salt concentration in mM.
    :type saltc: float.
    :param method: computation method to use. Only available method is
    'finnzymes'.
    :type method: str.

    '''

    # Universal gas constant
    R = 1.9872
    # Unit corrections for input paramters
    saltc = saltc / 1e3
    sc = 16.6 * log(saltc) / log(10.0)
    dnac = dnac / 1e9

    if method == 'finnzymes':
        params = finnzymes_par
    elif method == 'santalucia98':
        params = santalucia98_par
    else:
        '\'finnzymes\' is the only method that is currently supported'
    delta_H_par = params['delta_H']
    delta_S_par = params['delta_S']
    delta_H_par_err = params['delta_H_err']
    delta_S_par_err = params['delta_S_err']

    def collect_pairs(seq, pat):
    # This is faster than Biopython's overcount
        count = 0
        start = 0
        while True:
            start = seq.find(pat, start) + 1
            if start > 0:
                count += 1
            else:
                return count

    # shorter name
    seq = sequence
    seq = seq.upper()
    # Sum up the nearest-neighbor enthalpy and entropy counts*par
    H_keys = delta_H_par.keys()
    S_keys = delta_S_par.keys()
    dh = sum([collect_pairs(seq, x) * delta_H_par[x] for x in H_keys])
    ds = sum([collect_pairs(seq, x) * delta_S_par[x] for x in S_keys])

    # Error corrections
    # Initiation for seqs with a G-C pair vs only A-T
    if 'G' in seq or 'C' in seq:
        dh += delta_H_par_err['initGC']
        ds += delta_S_par_err['initGC']
    at_count = [x for x in seq if x == 'A' or x == 'T']
    if len(at_count) == len(seq):
        dh += delta_H_par_err['initAT']
        ds += delta_S_par_err['initAT']
    # 5' terminal T-A
    if seq.startswith('T'):
        dh += delta_H_par_err['5termT']
        ds += delta_S_par_err['5termT']
    # Correction for self-complementary sequences
    # The meaning of 'self-complementary' is not well-defined...
    if seq == reverse_complement(seq):
        dh += delta_H_par_err['symm']
        ds += delta_S_par_err['symm']

    # Unit corrections
    dh = dh * 1e3

    # These corrections are unaccounted for but are required
    # for the 'finnzymes' method to be accurate.
    if method == 'finnzymes':
        dh += 3400
        ds += 12.4
        tm = -dh / (R * log(dnac / 16) - ds) + sc - 273.15
    if method == 'santalucia98':
        k = dnac / 4.0
        ds = ds - 0.368 * (len(seq) - 1) * log(saltc)
        tm = -dh / ((R * log(k)) - ds) - 273.15

    return tm

finnzymes_par = {
    'delta_H': {
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
    'delta_H_err': {
        'initAT': 0.0,
        'initGC': 0.0,
        'symm': 0.0,
        '5termT': 0.0},
    'delta_S': {
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
    'delta_S_err': {
        'initAT': 0.0,
        'initGC': 0.0,
        'symm': 0.0,
        '5termT': 0.0}}

santalucia98_par = {
    'delta_H': {
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
    'delta_H_err': {
        'initAT': 0.0,
        'initGC': 0.0,
        'symm': 0.0,
        '5termT': -0.4},
    'delta_S': {
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
    'delta_S_err': {
        'initAT': 9.0,
        'initGC': 5.9,
        'symm': 1.4,
        '5termT': 0.0}}
