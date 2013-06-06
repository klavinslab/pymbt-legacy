'''
Module supplying methods/classes for generating overlapping oligo sequences
from a gene sequence.

'''

import csv
from math import floor
from pymbt.sequence_manipulation import reverse_complement
from pymbt.tm_calc import calc_tm
from pymbt.primer_design import design_primer_gene

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna


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
        self.seq = seq
        self.oligos = assembly_dict['oligos']
        self.overlaps = assembly_dict['overlaps']
        self.overlaps_tms = assembly_dict['overlaps_tms']
        self.overlaps_indices = assembly_dict['overlaps_indices']
        if primers:
            primers = design_primer_gene(seq, tm=primer_tm)
            self.primers = [x[0].upper() for x in primers]
            self.primer_tms = [x[1] for x in primers]

        # TODO: fix this problem automatically rather than warning
        for i in range(len(self.overlaps_indices) - 1):
            current_start = self.overlaps_indices[i + 1]
            current_end = self.overlaps_indices[i]
            if current_start <= current_end:
                print 'warning: overlapping overlaps!'

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
                o_tm = self.overlaps_tms[i]
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

    def write_map(self, path):
        starts = [self.seq.find(overlap) for overlap in
                  self.overlaps]
        overlap_lens = [len(overlap) for overlap in self.overlaps]
        features = []
        for i, start in enumerate(starts):
            location = FeatureLocation(ExactPosition(start),
                                       ExactPosition(start + overlap_lens[i]),
                                       strand=1)
            features.append(SeqFeature(location, type='misc_feature',
                            qualifiers={'label': ['overlap {}'.format(i)]}))
        seq_map = SeqRecord(Seq(self.seq, unambiguous_dna), features=features)
        SeqIO.write(seq_map, path, 'genbank')

    def __repr__(self):
        str1 = "An OligoAssembly consisting of "
        str2 = str(len(self.oligos)) + ' oligos.'
        return str1 + str2


def oligo_calc(seq, tm=72, length_range=(80, 200), require_even=True,
               start_5=True, oligo_number=None, overlap_min=20):
    '''Split a sequence into overlapping oligonucleotides with equal-Tm
    overlaps.

    :param seq: Input sequence (DNA).
    :type seq: str
    :param tm: Ideal Tm of the overlaps, in degrees C.
    :type tm: float
    :param length_range: Maximum oligo size (e.g. 60bp price point cutoff)
                         range - lower bound only matters if oligo_number
                         parameter is set.
    :type length_range: int 2-tuple
    :param require_even: Require that the number of oligonucleotides is even.
    :type require_even: bool
    :param start_5: Require that the first oligo's terminal side is 5\'.
    :type start_5: bool
    :param oligo_number: Attempt to build assembly with this many oligos,
                         starting with length_range min and incrementing by 10
                         up to length_range max.
    :type oligo_number: bool
    :param overlap_min: Minimum overlap size.
    :type overlap_min: int
    '''

    if len(seq) < length_range[0]:
        # If sequence can be built with just two oligos, do that
        oligos = [seq, reverse_complement(seq)]
        overlaps = [seq]
        overlaps_tms = [calc_tm(seq)]
        assembly_dict = {'oligos': oligos,
                         'overlaps': overlaps,
                         'overlaps_tms': overlaps_tms}
        return assembly_dict

    if oligo_number:
        # Make first attempt using length_range[0] and see what happens
        step = 3
        length_max = length_range[1]
        current_oligo_n = oligo_number + 1
        while current_oligo_n != oligo_number and length_max > length_range[0]:
            # Make oligos until target number is met. If impossible, make next
            # best thing
            # Tried starting with low range and going up - using too low
            # a 'max' results in bad things for longer sequences (overlaps
            # become longer than 80)

            # TODO: first run is redundant. fix it
            grown_overlaps = grow_overlaps(seq, tm, require_even, length_max,
                                           overlap_min)
            current_oligo_n = len(grown_overlaps[0])
            if current_oligo_n > oligo_number:
                break
            length_max -= step
    else:
        grown_overlaps = grow_overlaps(seq, tm, require_even, length_range[1],
                                       overlap_min)

    oligos, overlaps, overlaps_tms, overlaps_indices = grown_overlaps

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
                     'overlaps_tms': overlaps_tms,
                     'overlaps_indices': overlaps_indices}
    return assembly_dict


def grow_overlaps(seq, tm, require_even, length_max, overlap_min):
    '''
    :param seq: Input sequence (DNA).
    :type seq: str
    :param tm: Ideal Tm of the overlaps, in degrees C.
    :type tm: float
    :param require_even: Require that the number of oligonucleotides is even.
    :type require_even: bool
    :param length_max: Maximum oligo size (e.g. 60bp price point cutoff)
                       range.
    :type length_range: int
    :param overlap_min: Minimum overlap size.
    :type overlap_min: int

    '''
    # TODO: prevent growing overlaps from bumping into each other -
    # should halt when it happens, give warning, let user decide if they still
    # want the current construct
    # Another option would be to start over, moving the starting positions
    # near the problem region a little farther from each other - this would
    # put the AT-rich region in the middle of the spanning oligo

    # Number of oligos (initially) to start building
    oligo_n = int(floor(float(len(seq)) / length_max) + 1)

    # Adjust if even number of oligos is required
    if require_even:
        oligo_increment = 2
        if oligo_n % 2 == 1:
            oligo_n += 1
    else:
        oligo_increment = 1

    # Check to see if sequence can be spanned by that number of oligos given
    # a maximum length ceiling. If not, increase number of oligos.
    min_oligo_len = float(len(seq)) / oligo_n
    while min_oligo_len > length_max:
        oligo_n += oligo_increment
        min_oligo_len = float(len(seq)) / oligo_n

    # Initialize score list - list of tms. They start artificially low to make
    # min tm check simple.
    lowest_overlap_tm = -40000
    overlaps_tms = [lowest_overlap_tm] * (oligo_n - 1)

    # Loop until all overlaps meet minimum Tm
    tm_met = all([x >= tm for x in overlaps_tms])
    # TODO: len_unmet can loop forever or be exceedingly slow - e.g. generate
    # 900bp gene, tell it to do oligo_min of 50 and see what happens.
    len_met = False
    while(not tm_met or not len_met):
        # Calculate initial number of overlaps
        overlap_n = oligo_n - 1
        min_oligo_len = float(len(seq)) / oligo_n

        # Calculate overlap locations
        starts = [int(floor((x + 1) * min_oligo_len)) for x in
                  range(overlap_n)]
        ends = [index + 1 for index in starts]

        # Fencepost for while loop
        # Initial overlaps (1 base) and their tms
        overlaps = [seq[start:ends[i]] for i, start in enumerate(starts)]
        overlaps_tms = [calc_tm(overlap) for overlap in overlaps]
        lowest_overlap_tm = min(overlaps_tms)
        lowest_index = overlaps_tms.index(lowest_overlap_tm)
        # Initial oligos - includes the 1 base overlaps.
        # All the oligos are in the same direction - reverse
        # complementation of every other one happens later
        oligo_starts = [0] + starts
        oligo_ends = ends + [len(seq)]
        oligos = [seq[oligo_start:oligo_ends[i]] for i, oligo_start in
                  enumerate(oligo_starts)]

        # Initialize loop conditions
        maxed = [False for i in range(oligo_n)]

        while not (tm_met and len_met) and not any(maxed):
            # Recalculate overlaps and their Tms
            overlaps = [seq[oligo_starts[x + 1]:oligo_ends[x]] for x in
                        range(overlap_n)]
            overlaps_tms[lowest_index] = calc_tm(overlaps[lowest_index])

            # Find lowest-Tm overlap and its index.
            lowest_overlap_tm = min(overlaps_tms)
            lowest_index = overlaps_tms.index(lowest_overlap_tm)

            # Oligos to the left and right of the lowest Tm overlap
            oligo_to_left = len(oligos[lowest_index])
            oligo_to_right = len(oligos[lowest_index + 1])

            # If one of the oligos is max size, increase the other one
            if oligo_to_right == length_max:
                oligo_ends = _increase_right_oligo(oligo_ends, lowest_index)
            elif oligo_to_left == length_max:
                oligo_starts = _increase_left_oligo(oligo_starts, lowest_index)
            else:
                # Increase the smaller of the two oligos or,
                # if the same size, increase the one on the left
                # This biases the result but ensures it is deterministic
                if oligo_to_left > oligo_to_right:
                    oligo_starts = _increase_left_oligo(oligo_starts,
                                                        lowest_index)
                else:
                    oligo_ends = _increase_right_oligo(oligo_ends,
                                                       lowest_index)
            # Recalculate oligos from start and end indices
            oligos = [seq[oligo_starts[x]:oligo_ends[x]] for x in
                      range(oligo_n)]

            # Regenerate conditions
            maxed = [len(x) == length_max for x in oligos]
            tm_met = all([x >= tm for x in overlaps_tms])
            len_met = all([len(x) >= overlap_min for x in overlaps])

        oligo_n += oligo_increment
    overlaps_indices = [(oligo_starts[x + 1], oligo_ends[x]) for x in
                        range(overlap_n)]

    return oligos, overlaps, overlaps_tms, overlaps_indices


def _increase_right_oligo(positions_list, index):
    '''Increase oligo to the right.'''
    positions_list[index] += 1
    return positions_list


def _increase_left_oligo(positions_list, index):
    '''Increase oligo to the left.'''
    positions_list[index + 1] -= 1
    return positions_list
