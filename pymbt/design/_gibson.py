'''Gibson design module.'''
from pymbt import analysis, sequence
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
    # TODO: if sequence is too short (overlap len = seq len), raise exception
    dna1_primer = design_primer(dna1.reverse_complement(), **kwargs)
    dna2_primer = design_primer(dna2, **kwargs)
    if split == 'left':
        overhang_f = design_primer(dna1.reverse_complement(), tm=overlap_tm,
                                   tm_undershoot=0)
        overhang2 = overhang_f.primer().reverse_complement().flip()
        overhang1 = None
    elif split == 'right':
        overhang_r = design_primer(dna2, tm=overlap_tm, tm_undershoot=0)
        overhang1 = overhang_r.primer().reverse_complement().flip()
        overhang2 = None
    elif split == 'mixed':
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
        overhang1 = overlap_r.reverse_complement()
        overhang2 = overlap_l
    else:
        raise ValueError('split argument must be left, right, or mixed')
    primer1 = sequence.Primer(dna1_primer.anneal, tm=dna1_primer.tm,
                              overhang=overhang1)
    primer2 = sequence.Primer(dna2_primer.anneal, tm=dna2_primer.tm,
                              overhang=overhang2)

    return primer1, primer2


def gibson(seq_list, circular=True, splits='mixed', overlap_tm=65, **kwargs):
    '''Design Gibson primers given a set of sequences

    :param seq_list: List of DNA sequences to stitch together
    :type seq_list: list containing pymbt.sequence.DNA
    :param circular: If true, designs primers for making a circular construct.
                     If false, designs primers for a linear construct.
    :type circular: bool
    :param splits: Specifies locations of overlap. Must be either a single
                   entry of the same type as the 'split' parameter in
                   gibson_primers or a list of those types of the appropriate
                   length (for circular construct, len(seq_list), for
                   linear construct, len(seq_list) - 1)
    :type splits: str or list of str
    :param overlap_tm: Minimum Tm of overlap
    :type overlap_tm: float
    :param kwargs: keyword arguments to pass to design_primer
    :type kwargs: dict

    '''

    # Input checking
    if circular:
        n_overlaps = len(seq_list)
    else:
        n_overlaps = len(seq_list) - 1

    if type(splits) is str:
        splits = [splits] * n_overlaps
    else:
        if len(splits) != n_overlaps:
            raise ValueError("Incorrect number of 'splits' entries.")
        else:
            for split in splits:
                if split not in ['left', 'right', 'mixed']:
                    raise ValueError("Invalid 'splits' setting.")

    # If here, inputs were good
    # Design primers for linear constructs:
    primers_list = []
    for i, (left, right) in enumerate(zip(seq_list[:-1], seq_list[1:])):
        primers_list.append(gibson_primers(left, right, splits[i],
                                           overlap_tm=overlap_tm))
    if circular:
        primers_list.append(gibson_primers(seq_list[-1], seq_list[0],
                                           splits[-1], overlap_tm=overlap_tm))
    else:
        primer_f = design_primer(seq_list[0])
        primer_r = design_primer(seq_list[-1].reverse_complement())
        primers_list.append((primer_r, primer_f))

    # Primers are now in order of 'reverse for seq1, forward for seq2' config
    # Should be in 'forward and reverse primers for seq1, then seq2', etc
    # Just need to rotate one to the right
    flat = [y for x in primers_list for y in x]
    flat = [flat[-1]] + flat[:-1]
    grouped_primers = [(flat[2 * i], flat[2 * i + 1]) for i in
                       range(len(flat) / 2)]

    return grouped_primers
