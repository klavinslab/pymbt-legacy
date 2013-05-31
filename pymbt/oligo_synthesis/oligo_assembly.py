'''
Module supplying methods/classes for generating overlapping oligo sequences
from a gene sequence.

'''

import csv
from math import floor
from pymbt.sequence_manipulation import reverse_complement
from pymbt.tm_calc import calc_tm
from pymbt.primer_design import design_primer_gene


class OligoAssembly(object):
    '''Split a sequence into overlapping oligonucleotides with equal-Tm
    overlaps. Contains results and output methods.'''

    def __init__(self, seq, primers=False, primer_tm=60, **kwargs):
        '''
        :param seq: Input sequence (DNA).
        :type seq: str
        :param primers: Design cloning primers that bind the termini of the
                        input sequence.
        :type primers: bool
        :param primer_tm: Ideal Tm for the cloning primers, if applicable.
        :type primer_tm: float
        :param kwargs: Keyword arguments to pass to oligo_calc.

        '''

        assembly_dict = oligo_calc(seq=seq, **kwargs)
        self.oligos = assembly_dict['oligos']
        self.overlaps = assembly_dict['overlaps']
        self.overlap_tms = assembly_dict['overlap_tms']
        if primers:
            primers = design_primer_gene(seq, tm=primer_tm)
            self.primers = [x[0].upper() for x in primers]
            self.primer_tms = [x[1] for x in primers]

    def write(self, outpath):
        '''
        Write results out to a csv (comma-separated) file.

        :param outpath: path to csv file, including .csv extension.
        :type outpath: str

        '''

        oligo_writer = csv.writer(open(outpath, 'wb'), delimiter=',',
                                  quoting=csv.QUOTE_MINIMAL)
        oligo_writer.writerow(['name', 'oligo', 'notes'])
        for i, oligo in enumerate(self.oligos):
            name = 'oligo %i' % (i + 1)
            if i != len(self.oligos) - 1:
                o_tm = self.overlap_tms[i]
                len_tm = (len(oligo), o_tm)
                notes = 'oligo length: %d, overlap Tm: %.2f' % len_tm
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

    def __repr__(self):
        str1 = "An OligoAssembly consisting of "
        str2 = str(len(self.oligos)) + ' oligos.'
        return str1 + str2


def oligo_calc(seq, tm=72, length_range=(80, 200), require_even=True,
               start_5=True, oligo_number=None):
    '''Split a sequence into overlapping oligonucleotides with equal-Tm
    overlaps.

    :param seq: Input sequence (DNA).
    :type seq: str
    :param tm: Ideal Tm of the overlaps, in degrees C.
    :type tm: float
    :param length_range: Maximum oligo size (e.g. 60bp price point cutoff)
                         range.
    :type length_range: int 2-tuple
    :param require_even: Require that the number of oligonucleotides is even.
    :type require_even: bool
    :param start_5: Require that the first oligo's terminal side is 5\'.
    :type start_5: bool
    :param oligo_number: Attempt to build assembly with this many oligos,
                         starting with length_range min and incrementing by 10
                         up to length_range max.
    :type oligo_number: bool

    '''

    if len(seq) < length_range[0]:
        # If sequence can be built with just two oligos, do that
        oligos = [seq, reverse_complement(seq)]
        overlaps = [seq]
        overlap_tms = [calc_tm(seq)]
        assembly_dict = {'oligos': oligos,
                         'overlaps': overlaps,
                         'overlap_tms': overlap_tms}
        return assembly_dict

    if oligo_number:
        # Make first attempt using length_range[0] and see what happens
        step = 5
        length_max = length_range[0]
        current_oligo_n = oligo_number + 1
        while current_oligo_n != oligo_number and length_max < length_range[1]:
            # Make oligos until target number is met. If impossible, make next
            # best thing
            if current_oligo_n < oligo_number:
                break

            # TODO: first run is redundant. fix it
            grown_overlaps = grow_overlaps(seq, tm, require_even, length_max)
            current_oligo_n = len(grown_overlaps[0])
            length_max += step
    else:
        grown_overlaps = grow_overlaps(seq, tm, require_even, length_range[1])

    oligos, overlaps, overlap_tms = grown_overlaps

    if start_5:
        for i in [x for x in range(len(oligos)) if x % 2 == 1]:
            r_oligo = reverse_complement(oligos[i])
            oligos[i] = r_oligo
    else:
        for i in [x for x in range(len(oligos)) if x % 2 == 0]:
            r_oligo = reverse_complement(oligos[i])
            oligos[i] = r_oligo

    oligos = [x.upper() for x in oligos]
    assembly_dict = {'oligos': oligos,
                     'overlaps': overlaps,
                     'overlap_tms': overlap_tms}
    return assembly_dict


def grow_overlaps(seq, tm, require_even, length_max):
    oligo_n = int(floor(float(len(seq)) / length_max) + 1)

    if require_even:
        oligo_increment = 2
        if oligo_n % 2 == 1:
            oligo_n += 1
    else:
        oligo_increment = 1

    lowest_tm_overlap = -40000  # Arbitrary super negative number
    overlap_tms = [lowest_tm_overlap] * (oligo_n - 1)

    init = float(len(seq)) / oligo_n
    while init > length_max:
        oligo_n += oligo_increment
        init = float(len(seq)) / oligo_n

    # initial number of overlaps
    overlap_n = oligo_n - 1
    olap_iter = range(overlap_n)
    overlap_starts = [int(floor((x + 1) * init)) for x in olap_iter]
    overlap_ends = [x + 1 for x in overlap_starts]
    overlaps = [seq[overlap_starts[i]:overlap_ends[i]] for i in
                olap_iter]

    # Loop until all overlaps meet minimum Tm
    while(any([x <= tm for x in overlap_tms])):
        # Calculate:
        # initial number of overlaps
        # overlap locations
        overlap_n = oligo_n - 1
        olap_iter = range(overlap_n)
        init = float(len(seq)) / oligo_n
        overlap_starts = [int(floor((x + 1) * init)) for x in olap_iter]
        overlap_ends = [x + 1 for x in overlap_starts]
        overlaps = [seq[overlap_starts[i]:overlap_ends[i]] for i in
                    olap_iter]
        oligo_start = [0] + overlap_starts
        oligo_e = overlap_ends + [len(seq)]
        oligos = [seq[oligo_start[i]:oligo_e[i]] for i in olap_iter]

        # Get initial overlap tms (1 base?)
        overlap_tms = [calc_tm(x) for x in overlaps]

        # Initializations for next while loop

        # nonmaxed overlaps
        overlap_tms = [calc_tm(x) for x in overlaps]
        lowest_tm_overlap = min(overlap_tms)
        lowest_index = overlap_tms.index(lowest_tm_overlap)
        maxed = [False for i in range(oligo_n)]
        threshold_unmet = True

        while not all(maxed) and threshold_unmet and not any(maxed):
            overlaps = [seq[oligo_start[x + 1]:oligo_e[x]] for x in
                        olap_iter]
            # Calculate Tm of newly-increased overlap
            overlap_tms[lowest_index] = calc_tm(overlaps[lowest_index])
            # Identify overlap with the lowest Tm
            lowest_tm_overlap = min(overlap_tms)
            # Index of overlap with the lowest Tm
            lowest_index = overlap_tms.index(lowest_tm_overlap)
            # Oligos to the left and right of the lowest Tm overlap
            lowest_tm_overlap_left = len(oligos[lowest_index])
            lowest_tm_overlap_right = len(oligos[lowest_index + 1])

            def increase_right_oligo():
                '''Increase oligo to the right.'''
                oligo_e[lowest_index] += 1

            def increase_left_oligo():
                '''Increase oligo to the left.'''
                oligo_start[lowest_index + 1] -= 1

            # If one of the oligos is max size, increase the other one
            if lowest_tm_overlap_right == length_max:
                increase_right_oligo()
            elif lowest_tm_overlap_left == length_max:
                increase_left_oligo()
            else:
                # Increase the smaller of the two oligos or,
                # if the same size, increase the one on the left
                # This biases the result but ensures it is deterministic
                if lowest_tm_overlap_left > lowest_tm_overlap_right:
                    increase_left_oligo()
                else:
                    increase_right_oligo()

            # Recalculate oligos from start and end indices
            oligos = [seq[oligo_start[x]:oligo_e[x]] for x in
                      range(oligo_n)]
            maxed = [len(x) == length_max for x in oligos]
            threshold_unmet = any([x <= tm for x in overlap_tms])

        oligo_n += oligo_increment
    return oligos, overlaps, overlap_tms
