# Module supplying methods/classes for generating
# overlapping oligo sequences from a gene sequence.
# Input should be: a validated DNA sequence
# Output should: contain the set of oligos, the overlaps, and the overlap tms

#TODO:
# Clean up redundancies
# oligo_calc has terrible variable names
# may need to synchronize oligo calc class and function

import csv
from math import floor

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from pymbt.tm_calc import calc_tm
from pymbt.primer_design import design_primer_gene


class OligoAssembly(object):
    def __init__(self,
                 seq,
                 primers=False,
                 primer_tm=60,
                 **kwargs):

        assembly_dict = oligo_calc(seq=seq,
                                   **kwargs)
        self.oligos = assembly_dict['oligos']
        self.overlaps = assembly_dict['overlaps']
        self.overlap_tms = assembly_dict['olap_tms']
        if primers:
            primers = design_primer_gene(seq, tm=primer_tm)
            self.primers = [x[0].upper() for x in primers] 
            self.primer_tms = [x[1] for x in primers] 

    def __repr__(self):
        str1 = "An OligoAssembly consisting of "
        str2 = str(len(self.oligos)) + ' oligos.'
        return str1 + str2

    def write(self, outpath):
        oligo_writer = csv.writer(open(outpath, 'wb'),
                                  delimiter=',',
                                  quoting=csv.QUOTE_MINIMAL)
        oligo_writer.writerow(['name', 'oligo', 'notes'])
        for i, x in enumerate(self.oligos):
            name = 'oligo %i' % (i + 1)
            oligo = self.oligos[i]
            if i is not len(self.oligos) - 1:
                o_tm = self.overlap_tms[i]
                notes = 'oligo length: %d, overlap Tm: %.2f' % (len(oligo), o_tm)
            else:
                notes = 'oligo length: %d' % len(oligo)
            oligo_writer.writerow([name,
                                   oligo,
                                   notes])
        try:
            oligo_writer.writerow(['primer 1',
                                  self.primers[0],
                                  'Tm: %.2f' % self.primer_tms[0]])
            oligo_writer.writerow(['primer 2',
                                  self.primers[1],
                                  'Tm: %.2f' % self.primer_tms[1]])
        except AttributeError:
            pass
           


def oligo_calc(seq,
               tm=72,
               oligo_size=120,
               require_even=True,
               start_5=True,
               max_size=False,
               keep_even=True):

    if len(seq) < oligo_size:
        raise ValueError('Oligo size must be smaller than input sequence')
    # Ensures input type is valid
    seq = SeqRecord(Seq(seq, IUPAC.IUPACUnambiguousDNA)).seq.tostring()
    # Number of oligos to try first
    oligo_n = int(floor(float(len(seq)) / oligo_size) + 1)
    if require_even:
        oligo_increment = 2
        if oligo_n % 2 == 1:
            oligo_n += 1
    else:
        oligo_increment = 1
    lowest_tm = calc_tm('')
    olap_tms = [lowest_tm] * (oligo_n - 1)

    # Loop until all overlaps meet minimum Tm
    # (overlaps extended to full oligo_size)
    while(any([x <= tm for x in olap_tms])):
        # initial overlap locations
        init = float(len(seq)) / oligo_n
        o_l = [int(floor((x + 1) * init)) for x in range(oligo_n - 1)]
        o_start = [0] + o_l
        o_end = [x + 1 for x in o_l] + [len(seq)]

        # Initial data
        overlaps = [seq[o_start[i + 1]:o_end[i]] for i in range(oligo_n - 1)]
        olap_tms = [calc_tm(x) for x in overlaps]
        oligos = [seq[o_start[i]:o_end[i]] for i in range(oligo_n)]
        remain_olap = range(oligo_n - 1)
        remain_tm = [calc_tm(x) for x in overlaps]
        remaining_index = [i for i, x in enumerate(remain_tm)]
        lowest_tm = min(remain_tm)
        lowest_index = remain_tm.index(lowest_tm)
        maxed = [False for i in range(oligo_n)]
        threshold_unmet = True
        anymaxed = False

        while not all(maxed) and threshold_unmet and not anymaxed:
            o_iter = range(oligo_n - 1)
            overlaps = [seq[o_start[x + 1]:o_end[x]] for x in o_iter]
            olap_tms[lowest_index] = calc_tm(overlaps[lowest_index])
            remain_tm = [j for i, j in enumerate(olap_tms) if i in remain_olap]
            lowest_tm = min(remain_tm)
            lowest_index = remaining_index[remain_tm.index(lowest_tm)]
            relative_index = remain_tm.index(lowest_tm)
            oligo_l = len(oligos[lowest_index])
            oligo_r = len(oligos[lowest_index + 1])
            if oligo_l == oligo_size & oligo_r == oligo_size:
                remain_olap.pop(relative_index)
                remaining_index.pop(relative_index)
                changed = False
            elif oligo_r == oligo_size:
                changed = lowest_index
                o_end[changed] += 1
            elif oligo_l == oligo_size:
                changed = lowest_index + 1
                o_start[changed] -= 1
            else:
                # Increase smaller oligo or, if same size, the left one
                # Note: this also biases construction towards the left
                if oligo_l > oligo_r:
                    changed = lowest_index + 1
                    o_start[changed] -= 1
                else:
                    changed = lowest_index
                    o_end[changed] += 1

            oligos = [seq[o_start[x]:o_end[x]] for x in range(oligo_n)]
            maxed = [len(x) == oligo_size for x in oligos]

            if not max_size:
                threshold_unmet = any([x <= tm for x in olap_tms])

            if keep_even:
                anymaxed = any(maxed)

        oligo_n += oligo_increment

    if start_5:
        for i in [x for x in range(len(oligos)) if x % 2 == 1]:
            r_oligo = Seq(oligos[i]).reverse_complement().tostring()
            oligos[i] = r_oligo
    else:
        for i in [x for x in range(len(oligos)) if x % 2 == 0]:
            r_oligo = Seq(oligos[i]).reverse_complement().tostring()
            oligos[i] = r_oligo

    oligos = [x.upper() for x in oligos]
    assembly_dict = dict(oligos=oligos, overlaps=overlaps, olap_tms=olap_tms)
    return assembly_dict
