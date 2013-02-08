#TODO:
# sources: http://web.nmsu.edu/~snsm/classes/chem435/Lab4/
#   1. 'AT init' and 'GC init' terms - what is proper way to implement?
#   2. symmetry - factor it in? Most don't seem to.
#   3. Add the following methods - SantaLucia96, Sugimoto

'''Calculate the thermodynamic melting temperatures of nucleotide sequences
using the Finnzymes modified Breslauer 1986 parameters.'''

from math import log


def calc_tm(s, dnac=50, saltc=50, method='finnzymes'):
    '''Returns DNA/DNA tm using nearest neighbor thermodynamics.

    dnac is DNA concentration in nM
    saltc is salt concentration in mM.
    par is the parameter set to use (finnzymes, santaluciai, sugimoto)'''
    R = 1.9872  # universal gas constant.
    s = s.upper()

    def collect_pairs(seq, pat):
        count = 0
        start = 0
        while True:
            start = seq.find(pat, start) + 1
            if start > 0:
                count += 1
            else:
                return count

    if method == 'finnzymes':
        deltaH_par = finnzymes_par['deltaH']
        deltaS_par = finnzymes_par['deltaS']
    else:
        '\'finnzymes\' is the only currently supported method'

    # Sum up the nearest-neighbor enthalpy and entropy counts*par
    dh = sum([collect_pairs(s, x) * deltaH_par[x] for x in deltaH_par.keys()])
    ds = sum([collect_pairs(s, x) * deltaS_par[x] for x in deltaS_par.keys()])

    # unit corrections
    dh = dh * 1e3
    dh += 3400
    ds += 12.4
    saltc = saltc / 1e3
    dnac = dnac / 1e9

    # calculate salt contribution and Tm
    sc = 16.6 * log(saltc) / log(10.0)
    tm = -dh / (R * log(dnac / 16) - ds) + sc - 273.15

    return (tm)

finnzymes_par = {
    'deltaH': {
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
        'CC': 11.0,
        'termAT': 0,
        'termGC': 0
        },
    'deltaS': {
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
        'CC': 26.6,
        'termAT': 0,
        'termGC': 0
        },
    }

santalucia_par = {
    'deltaHParams': {
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
        'GG': 8,
        'CC': 8,
        'termAT': 2.3,
        'termGC': 0.1
        },
    'deltaSParams': {
        'AA': 22.2,
        'TT': 22.2,
        'AT': 20.4,
        'TA': 21.3,
        'CA': 22.7,
        'TG': 22.7,
        'GT': 22.4,
        'AC': 22.4,
        'CT': 21.0,
        'AG': 21.0,
        'GA': 22.2,
        'TC': 22.2,
        'CG': 27.2,
        'GC': 24.4,
        'GG': 19.9,
        'CC': 19.9,
        'termAT': 4.1,
        'termGC': -2.8
        }
    }
