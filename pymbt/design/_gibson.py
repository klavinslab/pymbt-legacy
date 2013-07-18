'''Gibson design module.'''
from pymbt import analysis
from pymbt.design import design_primer
# IDEA: Separate design of Gibson overlaps from design of primers


def gibson_primers(dna1, dna2, split, overlap_tm=65.0, **kwargs):
    '''Design Gibson primers given two DNA sequences (connect left to right)

    :param dna1: First piece of DNA for which to design primers. Once Gibsoned,
                 would be connected at its right side to dna2.
    :type dna1: pymbt.sequence.DNA
    :param dna2: First piece of DNA for which to design primers. Once Gibsoned,
                 would be connected at its right side to dna2.
    :type dna2: pymbt.sequence.DNA
    :param split: Specifies location of overlap. 'left' puts it on the 'dna1'
                  side (i.e. the primer to amplify dna2). 'right' puts it on
                  the dna2 side, and 'mixed' does a ~50:50 split
    :type split: str
    :param overlap_tm: Minimum Tm of overlap
    :type overlap_tm: float
    :param kwargs: keyword arguments to pass to design_primer
    :type kwargs: dict

    '''
    if split == 'left':
        overhang_f = design_primer(dna1.reverse_complement(), tm=overlap_tm,
                                   tm_undershoot=0)
        overhang = overhang_f.primer.reverse_complement()
        dna2_primer = design_primer(dna2, overhang=overhang, **kwargs)
        dna1_primer = design_primer(dna1.reverse_complement(), **kwargs)
    elif split == 'right':
        overhang_r = design_primer(dna2, tm=overlap_tm, tm_undershoot=0)
        overhang = overhang_r.primer.reverse_complement()
        dna2_primer = design_primer(dna2, **kwargs)
        dna1_primer = design_primer(dna1.reverse_complement(),
                                    overhang=overhang,
                                    **kwargs)
    elif split == 'mixed':
        # TODO: the overlaps are wrong!
        overlap_l = dna1[0:0]
        overlap_r = dna2[0]
        overlap_melt = analysis.tm(overlap_r)
        while overlap_melt < overlap_tm:
            rlen = len(overlap_r)
            llen = len(overlap_l)
            if rlen > llen:
                overlap_l = dna1[-(rlen + 1):]
            else:
                overlap_r = dna2[:(llen + 1)]
            overlap = overlap_l + overlap_r
            overlap_melt = analysis.tm(overlap)
        dna2_primer = design_primer(dna2, overhang=overlap_l, **kwargs)
        dna1_primer = design_primer(dna1.reverse_complement(),
                                    overhang=overlap_r.reverse_complement(),
                                    **kwargs)
    else:
        raise ValueError('split argument must be left, right, or mixed')

    return dna1_primer, dna2_primer
