'''Primer design tools.'''

from pymbt.tm_calc import calc_tm
from pymbt.sequence_manipulation import check_alphabet
from pymbt.sequence_manipulation import reverse_complement


def design_primer(seq, tm=72, min_len=10, tm_undershoot=1, tm_overshoot=3,
                  end_gc=False, tm_method='finnzymes', overhang=''):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

    :param seq: Input sequence.
    :type seq: str
    :param tm: Ideal primer Tm in degrees C.
    :type tm: float
    :param min_len: Minimum primer length.
    :type min_len: int
    :param tm_undershoot: Allowed Tm undershoot.
    :type tm_undershoot: float
    :param tm_overshoot: Allowed Tm overshoot.
    :type tm_overshoot: float
    :param end_gc: Obey the 'end on G or C' rule.
    :type end_gc: bool
    :param tm_method: Melting temp calculator method to use.
    :type tm_method: string
    :param overhang: Append the primer to this overhang sequence.
    :type overhang: str

    '''

    # Check Tm of input sequence to see if it's already too low
    seq_tm = calc_tm(seq, method=tm_method)
    if seq_tm < tm - tm_undershoot:
        err = 'Input sequence Tm is lower than primer Tm setting'
        raise Exception(err)

    max_tm = tm + tm_overshoot
    bases = min_len
    # Trim down max length to increase efficiency
    # Pretty much impossible for an annealing sequence to need more than 90bp
    seq = seq[0:90]

    # First, generate a range of primers and Tms:
    #     range: from min_len to 'tm' + tm_overshoot
    primers = []
    tms = []
    primer_tm = 0
    primer_len = 0
    while primer_tm <= max_tm and primer_len <= 80 and primer_len != len(seq):
        new_primer = seq[0:bases]
        new_tm = calc_tm(new_primer, method=tm_method)
        primers.append(new_primer)
        tms.append(new_tm)
        bases += 1
        primer_tm = tms[-1]
        primer_len = len(new_primer)

    # Trim primer list based on tm_undershoot and end_gc
    terr = tm - tm_undershoot
    #print terr
    #print tms
    primers = [primers[i] for i, x in enumerate(tms) if x >= terr]
    tms = [x for i, x in enumerate(tms) if x >= terr]

    if end_gc:
        primers = [x for x in primers if x.endswith(('C', 'G'))]
        tms = [tms[i] for i, x in enumerate(primers) if x.endswith(('C', 'G'))]

    if not primers:
        raise Exception('No primers could be generated using these settings')

    # Find the primer closest to the set Tm
    tm_diffs = [abs(x - tm) for x in tms]
    best_index = tm_diffs.index(min(tm_diffs))
    best_primer = (primers[best_index], tms[best_index])

    if overhang:
        check_alphabet(overhang)
    best_primer = (overhang + best_primer[0], best_primer[1])

    return best_primer


def design_primer_gene(seq, overhangs=None, **kwargs):
    '''
    Design both forward and reverse primers for a sequence.

    :param seq: Input sequence.
    :type seq: str
    :param overhangs: List of overhang sequences.
    :type overhangs: list
    :param kwargs: Keyword arguments to feed to design_primer.

    '''
    if not overhangs:
        overhangs = ['', '']

    fwd_primer = design_primer(seq,
                               overhang=overhangs[0],
                               **kwargs)
    rev_primer = design_primer(reverse_complement(seq),
                               overhang=overhangs[1],
                               **kwargs)
    return [fwd_primer, rev_primer]
