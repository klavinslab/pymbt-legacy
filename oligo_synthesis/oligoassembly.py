# Module supplying methods/classes for generating
# overlapping oligo sequences from a gene sequence.
# Input should be: a SeqRecord (this ensures it's a DNA sequence)
# Output should: contain the set of oligos, the overlaps, and the overlap tms
# TODO:
# may need to synchronize oligo calc class and function

import csv
from math import ceil

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from pymbt.calc_tm import calc_tm


class OligoAssembly(object):
    # Primary class, return all useful data (and oligos)
    def __init__(self,
                 seq,
                 tm=65,
                 oligo_size=120,
                 require_even=True,
                 start_5=True):

        results = oligo_calc(seq, tm, oligo_size, require_even, start_5)

        self.oligos = results['oligos']
        self.olaps = results['overlaps']
        self.tms = results['overlap_tms']

    def __repr__(self):
        str1 = "An OligoAssembly consisting of "
        str2 = str(len(self.oligos)) + ' oligos.'
        return(str1 + str2)

    def write(self, outpath):
        oligo_writer = csv.writer(open(outpath, 'wb'),
                                  delimiter=',',
                                  quoting=csv.QUOTE_NONNUMERIC)
        oligo_writer.writerow(['oligo', 'overlap', 'overlap_tm'])
        for i in range(len(self.oligos) - 1):
            oligo_writer.writerow([self.oligos[i], self.olaps[i], self.tms[i]])
        oligo_writer.writerow([self.oligos[len(self.oligos) - 1], 'NA', 'NA'])


def oligo_calc(seq, tm=65, oligo_size=120, require_even=False, start_5=True):
    # Setup - handling input and doing initial calculations
    if len(seq) < oligo_size:
        raise ValueError('Oligo size must be smaller than input sequence')
    seq = SeqRecord(Seq(seq, IUPAC.IUPACUnambiguousDNA))
    seq = seq.seq.tostring()
    #Absolute minimum oligos
    oligo_n = int(ceil(float(len(seq)) / oligo_size))
    if require_even:
        oligo_increment = 2
        if oligo_n % 2 == 1:
            oligo_n += 1
    else:
        oligo_increment = 1

    min_tm = calc_tm('')
    overlap_tms = [min_tm] * (oligo_n - 1)

    # Loop until all overlaps meet minimum Tm
    # (overlaps extended to full oligo_size)
    while(any([x <= tm for x in overlap_tms])):
        overlap_tms = [min_tm] * (oligo_n - 1)
        print('Trying with ' + str(oligo_n) + ' oligos.')
        # initial overlap locations
        o_range = range(oligo_n - 1)
        o_l = [int((x + 1) * float(len(seq)) / (oligo_n)) for x in o_range]
        # overlap starts
        o_s = [0] + o_l
        # overlap ends
        o_e = [x + 1 for x in o_l] + [len(seq)]
        # calculate initial overlaps
        overlaps = [seq[o_s[x + 1]:o_e[x]] for x in range(oligo_n - 1)]
        # Initial set of oligos
        oligos = [seq[o_s[x]:o_e[x]] for x in range(oligo_n)]
        # No oligos are maxed at first
        maxed = [False for x in range(oligo_n)]
        # No overlaps are maxed at first
        o_nomax = range(oligo_n - 1)
        while not all(maxed):
            # Find lowest overlap Tm index location
            # If more than one entry, it just picks the earliest one
            # Ignores any that are maxed
            nomax_tms = [i for j, i in enumerate(overlap_tms) if j in o_nomax]
            min_tm = overlap_tms.index(min(nomax_tms))
            oligo_l = len(oligos[min_tm])
            oligo_r = len(oligos[min_tm + 1])

            # Now decide which oligo (left or right) to increase
            if oligo_l == oligo_size & oligo_r == oligo_size:
                # If both have met the threshold, ignore that overlap
                o_nomax.pop(o_nomax.index(min_tm))
            elif oligo_r == oligo_size:
                # If right oligo is already maxed out, increase left
                changed = min_tm
                o_e[changed] += 1
            elif oligo_l == oligo_size:
                # If left oligo is already maxed out, increase right
                changed = min_tm + 1
                o_s[changed] -= 1
            else:
                # Increase smaller oligo or, if same size, the first one
                if oligo_l > oligo_r:
                    changed = min_tm + 1
                    o_s[changed] -= 1
                else:
                    changed = min_tm
                    o_e[changed] += 1

            # Recalculate oligos
            oligos = [seq[o_s[x]:o_e[x]] for x in range(oligo_n)]
            # Recalculate maxed-out oligos
            maxed = [len(x) == oligo_size for x in oligos]
            # Recalculate all overlaps
            overlaps = [seq[o_s[x + 1]:o_e[x]] for x in range(oligo_n - 1)]
            # Find the Tm of the changed oligo
            overlap_tms[min_tm] = calc_tm(overlaps[min_tm])

        oligo_n += oligo_increment

    # Reverse complement every second oligo
    for i in [x for x in range(len(oligos)) if x % 2 == 1]:
        r_oligo = Seq(oligos[i]).reverse_complement().tostring()
        oligos[i] = r_oligo
    result = dict(oligos=oligos, overlaps=overlaps, overlap_tms=overlap_tms)
    return(result)
