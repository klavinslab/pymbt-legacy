from pymbt.tm_calc import calc_tm
from pymbt.sequence_manipulation import check_alphabet
from pymbt.sequence_manipulation import reverse_complement


def design_primer(seq,
                  tm=72,
                  minstart=10,
                  tm_errorminus=1,
                  tm_errorplus=3,
                  endGC=False,
                  tail=''):
    # Check Tm of input sequence to see if it's already too low
    seq_tm = calc_tm(seq)
    if seq_tm < tm - tm_errorminus:
        err = 'Input sequence Tm is lower than primer Tm setting'
        raise Exception(err)

    max_tm = tm + tm_errorplus
    min_tm = tm - tm_errorminus
    bases = minstart
    # Trim down max length to increase efficiency
    # Pretty much impossible for an annealing sequence to need more than 90bp
    seq = seq[0:90]

    # First, generate a range of primers and Tms:
    #     range: from minstart to 'tm' + tm_errorplus
    primers = []
    tms = []
    primer_tm = 0
    while primer_tm <= max_tm:
        new_primer = seq[0:bases]
        new_tm = calc_tm(new_primer)
        primers.append(new_primer)
        tms.append(new_tm)
        bases += 1
        primer_tm = tms[-1]

    # Trim primer list based on tm_errorminus and endGC
    terr = tm - tm_errorminus
    primers = [primers[i] for i, x in enumerate(tms) if x >= terr]
    tms = [x for i, x in enumerate(tms) if x >= terr]
    if endGC:
        primers = [x for x in primers if x.endswith(('C', 'G'))]
        tms = [tms[i] for i, x in enumerate(primers) if x.endswith(('C', 'G'))]
    if not primers:
        raise Exception('No primers could be generated using these settings')

    # Find the primer closest to the set Tm
    tm_diffs = [abs(x - tm) for x in tms]
    best_index = tm_diffs.index(min(tm_diffs))
    best_primer = (primers[best_index], tms[best_index])

    if tail:
        check_alphabet(tail)
    best_primer = (tail + best_primer[0], best_primer[1])

    return best_primer


def design_primer_gene(seq,
                       tails=['', ''],
                       **kwargs):
    fwd_primer = design_primer(seq,
                               tail=tails[0],
                               **kwargs)
    rev_primer = design_primer(reverse_complement(seq),
                               tail=tails[1],
                               **kwargs)
    return [fwd_primer, rev_primer]
