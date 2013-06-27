'''
Module supplying methods/classes for generating overlapping oligo sequences
from a gene sequence.

'''

import csv
from math import floor
from pymbt import analysis
from pymbt.design import DesignPrimerGene

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna


class OligoAssembly(object):
    '''Split a sequence into overlapping oligonucleotides with equal-Tm
    overlaps. Contains results and output methods.'''

    def __init__(self, seq, tm=72, length_range=(80, 200), require_even=True,
                 start_5=True, oligo_number=None, overlap_min=20,
                 min_exception=False, primers=False, primer_tm=60):

        '''
        :param seq: Input sequence (DNA).
        :type seq: DNA object
        :param primers: Design cloning primers that bind the termini of the
                        input sequence.
        :type primers: bool
        :param primer_tm: Ideal Tm for the cloning primers, if applicable.
        :type primer_tm: float
        :param length_range: Maximum oligo size (e.g. 60bp price point cutoff)
                             range - lower bound only matters if oligo_number
                             parameter is set.
        :type length_range: int 2-tuple
        :param require_even: Require that the number of oligonucleotides is
                             even.
        :type require_even: bool
        :param start_5: Require that the first oligo's terminal side is 5\'.
        :type start_5: bool
        :param oligo_number: Attempt to build assembly with this many oligos,
                             starting with length_range min and incrementing by
                             10 up to length_range max.
        :type oligo_number: bool
        :param overlap_min: Minimum overlap size.
        :type overlap_min: int
        :param min_exception: In order to meet tm and overlap_min requirements,
                              allow overlaps less than overlap_min to continue
                              growing above tm setpoint.
        :type min_exception: bool

        '''

        self.kwargs = {'tm': tm, 'length_range': length_range,
                       'require_even': require_even, 'start_5': start_5,
                       'oligo_number': oligo_number,
                       'overlap_min': overlap_min,
                       'min_exception': min_exception}

        self.template = seq

        self.settings = {'primers': True, 'primer_tm': primer_tm}
        self._has_run = False

    def run(self):
        '''
        Execute design.

        '''

        # Input parameters needed to design the oligos
        length_range = self.kwargs['length_range']
        oligo_number = self.kwargs['oligo_number']
        require_even = self.kwargs['require_even']
        tm = self.kwargs['tm']
        overlap_min = self.kwargs['overlap_min']
        min_exception = self.kwargs['min_exception']
        start_5 = self.kwargs['start_5']

        if len(self.template) < length_range[0]:
            # If sequence can be built with just two oligos, do that
            oligos = [self.template, self.template.reverse_complement()]
            overlaps = [self.template]
            overlap_tms = [analysis.Tm(self.template).run()]
            assembly_dict = {'oligos': oligos, 'overlaps': overlaps,
                             'overlap_tms': overlap_tms}

            return assembly_dict

        if oligo_number:
            # Make first attempt using length_range[1] and see what happens
            step = 3  # Decrease max range by this amount each iteration

            length_max = length_range[1]
            current_oligo_n = oligo_number + 1
            oligo_n_met = False
            above_min_len = length_max > length_range[0]
            while not oligo_n_met and above_min_len:
                # Starting with low range and going up doesnt work for longer
                # sequence (overlaps become longer than 80)
                assembly = _grow_overlaps(self.template, tm, require_even,
                                          length_max, overlap_min,
                                          min_exception)
                current_oligo_n = len(assembly[0])
                if current_oligo_n > oligo_number:
                    break
                length_max -= step
                oligo_n_met = current_oligo_n == oligo_number
        else:
            assembly = _grow_overlaps(self.template, tm, require_even,
                                      length_range[1], overlap_min,
                                      min_exception)

        oligos, overlaps, overlap_tms, overlap_indices = assembly

        if start_5:
            for i in [x for x in range(len(oligos)) if x % 2 == 1]:
                oligos[i] = oligos[i].reverse_complement()
        else:
            for i in [x for x in range(len(oligos)) if x % 2 == 0]:
                oligos[i] = oligos[i].reverse_complement()

        # Make single-stranded
        oligos = [oligo.set_stranded('ss') for oligo in oligos]

        assembly_dict = {'oligos': oligos,
                         'overlaps': overlaps,
                         'overlap_tms': overlap_tms,
                         'overlap_indices': overlap_indices}

        self.oligos = assembly_dict['oligos']
        self.overlaps = assembly_dict['overlaps']
        self.overlap_tms = assembly_dict['overlap_tms']
        self.overlap_indices = assembly_dict['overlap_indices']

        if self.settings['primers']:
            primer_design = DesignPrimerGene(self.template,
                                             tm=self.settings['primer_tm'])
            self.primers_tms = primer_design.run()
            self.primers = [x[0] for x in self.primers_tms]
            self.primer_tms = [x[1] for x in self.primers_tms]

        for i in range(len(self.overlap_indices) - 1):
            current_start = self.overlap_indices[i + 1][0]
            current_end = self.overlap_indices[i][1]
            if current_start <= current_end:
                print 'warning: overlapping overlaps!'

        self._has_run = True

    def write(self, outpath):
        '''
        Write results out to a csv (comma-separated) file.

        :param outpath: path to csv file, including .csv extension.
        :type outpath: str

        '''

        with open(outpath, 'wb') as oligo_file:
            oligo_writer = csv.writer(oligo_file, delimiter=',',
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
        if self._has_run:
            str1 = "An OligoAssembly consisting of "
            str2 = str(len(self.oligos)) + ' oligos.'
            return str1 + str2
        else:
            return 'An OligoAssembly instance that has not been run yet.'


def _grow_overlaps(seq, tm, require_even, length_max, overlap_min,
                   min_exception):
    '''
    :param seq: Input sequence (DNA).
    :type seq: DNA object
    :param tm: Ideal Tm of the overlaps, in degrees C.
    :type tm: float
    :param require_even: Require that the number of oligonucleotides is even.
    :type require_even: bool
    :param length_max: Maximum oligo size (e.g. 60bp price point cutoff)
                       range.
    :type length_range: int
    :param overlap_min: Minimum overlap size.
    :type overlap_min: int
    :param min_exception: In order to meet tm and overlap_min requirements,
                          allow overlaps less than overlap_min to continue
                          growing above tm setpoint.
    :type min_exception: bool

    '''

    # TODO: prevent growing overlaps from bumping into each other -
    # should halt when it happens, give warning, let user decide if they still
    # want the current construct
    # Another option would be to start over, moving the starting positions
    # near the problem region a little farther from each other - this would
    # put the AT-rich region in the middle of the spanning oligo

    # Try bare minimum number of oligos
    oligo_n = int(floor(float(len(seq)) / length_max) + 1)

    # Adjust number of oligos if even number required
    if require_even:
        oligo_increment = 2
        if oligo_n % 2 == 1:
            oligo_n += 1
    else:
        oligo_increment = 1

    # Double check whether current oligos can even span the sequence
    oligo_len = float(len(seq)) / oligo_n
    while oligo_len > length_max:
        oligo_n += oligo_increment
        oligo_len = float(len(seq)) / oligo_n

    # Loop until all overlaps meet minimum Tm and length
    tm_met = False
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
        oligo_index = [oligo_starts, oligo_ends]

        oligos = [seq[start:end] for start, end in zip(*oligo_index)]

        # Oligo won't be maxed in first pass. tm_met and len_met will be false
        maxed = False

        while not (tm_met and len_met) and not maxed:
            # Recalculate overlaps and their Tms
            overlaps = _recalculate_overlaps(seq, overlaps, oligo_index)
            # Tm calculation is bottleneck - only recalculate changed overlap
            overlap_tms[index] = analysis.Tm(overlaps[index]).run()
            # Find lowest-Tm overlap and its index.
            index = overlap_tms.index(min(overlap_tms))
            # Move overlap at that index
            oligos = _expand_overlap(seq, oligo_index, index, oligos,
                                     length_max)
            # Regenerate conditions
            maxed = any([len(x) == length_max for x in oligos])
            tm_met = all([x >= tm for x in overlap_tms])
            if min_exception:
                len_met = True
            else:
                len_met = all([len(x) >= overlap_min for x in overlaps])

        # TODO: add test for min_exception case (use rob's sequence from
        # 20130624 with 65C Tm)
        if min_exception:
            len_met = all([len(x) >= overlap_min for x in overlaps])

            # See if len_met is true - if so do nothing
            if len_met:
                break
            else:
                while not len_met and not maxed:
                    # Recalculate overlaps and their Tms
                    overlaps = _recalculate_overlaps(seq, overlaps,
                                                     oligo_index)
                    # Overlap to increase is the shortest one
                    overlap_lens = [len(overlap) for overlap in overlaps]
                    index = overlap_lens.index(min(overlap_lens))
                    # Increase left or right oligo
                    oligos = _expand_overlap(seq, oligo_index, index, oligos,
                                             length_max)
                    # Recalculate conditions
                    maxed = any([len(x) == length_max for x in oligos])
                    len_met = all([len(x) >= overlap_min for x in overlaps])

                # Recalculate tms to reflect any changes (some are redundant)
                overlap_tms[index] = analysis.Tm(overlaps[index]).run()

                # Outcome could be that len_met happened *or* maxed out
                # length of one of the oligos. If len_met happened, should be
                # done so long as tm_met has been satisfied. If maxed happened,
                # len_met will not have been met, even if tm_met is satisfied,
                # and script will reattempt with more oligos

        oligo_n += oligo_increment

    # Calculate location of overlaps
    overlap_indices = [(oligo_index[0][x + 1], oligo_index[1][x]) for x in
                       range(overlap_n)]

    return oligos, overlaps, overlap_tms, overlap_indices


def _recalculate_overlaps(seq, overlaps, oligo_index):
    for i, overlap in enumerate(overlaps):
        first_index = oligo_index[0][i + 1]
        second_index = oligo_index[1][i]
        overlap = seq[first_index:second_index]
        overlaps[i] = overlap

    return overlaps


def _expand_overlap(seq, oligo_index, index, oligos, length_max):
    '''
    Given an overlap to increase, in which direction should it increase?

    :param oligo_index: index of oligo starts and stops
    :type oligo_index: list
    :param index: index of the oligo
    :type index: int
    :param left_len: length of left oligo
    :type left_len: int
    :param right_len: length of right oligo
    :type right_len: int
    :param length_max: length ceiling
    :type length_max: int

    '''

    left_len = len(oligos[index])
    right_len = len(oligos[index + 1])

    # If one of the oligos is max size, increase the other one
    if right_len == length_max:
        oligo_index[1] = _adjust_overlap(oligo_index[1], index, 'right')
    elif left_len == length_max:
        oligo_index[0] = _adjust_overlap(oligo_index[0], index, 'left')
    else:
        if left_len > right_len:
            oligo_index[0] = _adjust_overlap(oligo_index[0], index, 'left')
        else:
            oligo_index[1] = _adjust_overlap(oligo_index[1], index, 'right')

    # Recalculate oligos from start and end indices
    oligos = [seq[start:end] for start, end in zip(*oligo_index)]

    return oligos


def _adjust_overlap(positions_list, index, direction):
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
