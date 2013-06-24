'''
Module supplying methods/classes for generating overlapping oligo sequences
from a gene sequence.

'''

import csv
from math import floor
from pymbt.sequence.utils import reverse_complement
from pymbt import analysis
from pymbt.design import DesignPrimerGene
from pymbt.sequence.utils import check_instance

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna

# TODO: note - this is waaaaaaaay slower using DNA cloning language than
# using raw strings. It's spending tons of time checking the alphabet
# of something


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

        self.template = seq
        check_instance(self.template)
        assembly_dict = oligo_calc(seq=self.template, **kwargs)

        self.oligos = assembly_dict['oligos']
        self.overlaps = assembly_dict['overlaps']
        self.overlap_tms = assembly_dict['overlap_tms']
        self.overlap_indices = assembly_dict['overlap_indices']
        if primers:
            primers = DesignPrimerGene(seq, tm=primer_tm).run()
            self.primers = [x[0].upper() for x in primers]
            self.primer_tms = [x[1] for x in primers]

        # TODO: fix this problem automatically rather than warning
        for i in range(len(self.overlap_indices) - 1):
            current_start = self.overlap_indices[i + 1][0]
            current_end = self.overlap_indices[i][1]
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
            name = 'oligo {}'.format(i + 1)
            oligo_len = len(oligo)
            if i != len(self.oligos) - 1:
                oligo_tm = self.overlap_tms[i]
                notes = 'oligo length: {}, '.format(oligo_len) + \
                        'overlap Tm: {:.2f}'.format(oligo_tm)
            else:
                notes = 'oligo length: {}'.format(oligo_len)
            oligo_writer.writerow([name,
                                   oligo,
                                   notes])
        try:
            oligo_writer.writerow(['primer 1',
                                  self.primers[0],
                                  'Tm: {:.2f}'.format(self.primer_tms[0])])
            oligo_writer.writerow(['primer 2',
                                  self.primers[1],
                                  'Tm: {:.2f}'.format(self.primer_tms[1])])
        except AttributeError:
            pass

    def write_map(self, path):
        '''
        Write genbank map highlighting overlaps to file.

        :param path: full path to .gb file to write.
        :type path: str

        '''
        starts = [index[0] for index in self.overlap_indices]
        overlap_lens = [len(overlap) for overlap in self.overlaps]
        features = []
        for i, start in enumerate(starts):
            location = FeatureLocation(ExactPosition(start),
                                       ExactPosition(start + overlap_lens[i]),
                                       strand=1)
            overlap_name = 'overlap {}'.format(i + 1)
            features.append(SeqFeature(location, type='misc_feature',
                            qualifiers={'label': [overlap_name]}))
        seq_map = SeqRecord(Seq(self.template, unambiguous_dna),
                            features=features)
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
        oligos = [seq, reverse_complement(seq, 'dna')]
        overlaps = [seq]
        overlap_tms = [analysis.Tm(seq).run()]
        assembly_dict = {'oligos': oligos,
                         'overlaps': overlaps,
                         'overlap_tms': overlap_tms}
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

    oligos, overlaps, overlap_tms, overlap_indices = grown_overlaps

    if start_5:
        for i in [x for x in range(len(oligos)) if x % 2 == 1]:
            r_oligo = oligos[i].reverse_complement()
            #reverse_complement(oligos[i])
            oligos[i] = r_oligo
    else:
        for i in [x for x in range(len(oligos)) if x % 2 == 0]:
            r_oligo = oligos[i].reverse_complement()
            oligos[i] = r_oligo

    oligos = [str(x) for x in oligos]
    assembly_dict = {'oligos': oligos,
                     'overlaps': overlaps,
                     'overlap_tms': overlap_tms,
                     'overlap_indices': overlap_indices}
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

    #
    # Number of oligos (initially) to start building
    #

    # Bare minimum - usually not enough oligos
    oligo_n = int(floor(float(len(seq)) / length_max) + 1)

    # Adjust if even number of oligos is required
    if require_even:
        oligo_increment = 2
        if oligo_n % 2 == 1:
            oligo_n += 1
    else:
        oligo_increment = 1

    # Check to see if sequence can be spanned by that number of oligos given a
    # maximum length ceiling. If not, increase number of oligos.
    oligo_len = float(len(seq)) / oligo_n
    while oligo_len > length_max:
        oligo_n += oligo_increment
        oligo_len = float(len(seq)) / oligo_n

    # Initialize score list - list of tms. They start artificially low to make
    # min tm check simple.
    overlap_tms = [-40000] * (oligo_n - 1)

    # Loop until all overlaps meet minimum Tm
    tm_met = all([x >= tm for x in overlap_tms])

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
        overlaps = [seq[start:end] for start, end in zip(starts, ends)]
        overlap_tms = [analysis.Tm(overlap).run() for overlap in overlaps]
        index = overlap_tms.index(min(overlap_tms))
        # Initial oligos - includes the 1 base overlaps.
        # All the oligos are in the same direction - reverse
        # complementation of every other one happens later
        oligo_starts = [0] + starts
        oligo_ends = ends + [len(seq)]
        oligos = [seq[start:end] for start, end in
                  zip(oligo_starts, oligo_ends)]

        maxed = False

        while not (tm_met and len_met) and not maxed:
            # Recalculate overlaps and their Tms
            overlaps = [seq[oligo_starts[x + 1]:oligo_ends[x]] for x in
                        range(overlap_n)]
            overlap_tms[index] = analysis.Tm(overlaps[index]).run()

            # Find lowest-Tm overlap and its index.
            index = overlap_tms.index(min(overlap_tms))

            # Oligos to the left and right of the lowest Tm overlap
            left_oligo = len(oligos[index])
            right_oligo = len(oligos[index + 1])

            # If one of the oligos is max size, increase the other one
            if right_oligo == length_max:
                oligo_ends = _increase_overlap(oligo_ends, index, 'right')
            elif left_oligo == length_max:
                oligo_starts = _increase_overlap(oligo_starts, index, 'left')
            else:
                if left_oligo > right_oligo:
                    oligo_starts = _increase_overlap(oligo_starts, index,
                                                     'left')
                else:
                    oligo_ends = _increase_overlap(oligo_ends, index, 'right')

            # Recalculate oligos from start and end indices
            oligos = [seq[start:end] for start, end in
                      zip(oligo_starts, oligo_ends)]

            # Regenerate conditions
            maxed = any([len(x) == length_max for x in oligos])
            tm_met = all([x >= tm for x in overlap_tms])
            len_met = all([len(x) >= overlap_min for x in overlaps])

        oligo_n += oligo_increment

    overlap_indices = [(oligo_starts[x + 1], oligo_ends[x]) for x in
                       range(overlap_n)]

    return oligos, overlaps, overlap_tms, overlap_indices


def _increase_overlap(positions_list, index, direction):
    '''
    Increase overlap to the right or left of an index.

    :param positions_list: list of overlap positions
    :type positions_list: list
    :param index: index of the overlap to increase.
    :type index: int
    :param direction: which side of the overlap to increase - left or right.
    :type direction: str

    '''

    if direction == 'left':
        positions_list[index + 1] -= 1
    elif direction == 'right':
        positions_list[index] += 1
    else:
        raise ValueError('direction must be \'left\' or \'right\'.')

    return positions_list
