'''Sanger sequencing alignment tools.'''

import os
import math
from matplotlib import pyplot
from matplotlib import cm
from Bio import SeqIO
from pymbt.sequence_analysis.needle import needle
from pymbt.sequence_manipulation import reverse_complement

# TODO:
# consensus / master sequence for plotting / report / analysis
# too much work is done initializing Sanger. Move calculations to methods.


class Sanger(object):
    '''Align and analyze Sanger sequencing results.'''
    def __init__(self, ref, res):
        '''
        :param ref: Reference sequence.
        :type ref: str
        :param res: Sequencing result string. A list of strings is also valid.
        :type res: str

        '''

        # make results a list if there's just one
        if type(res) == str:
            res = [res]
        self.ref_raw = ref
        self.resnames = [x.name for x in res]
        #self.resnames = [x.description for x in res]
        ref = ref.seq.tostring()
        res = [x.seq.tostring() for x in res]
        # Reduce results to largest unambiguous segment and align
        split = [x.split('N') for x in res]
        lengths = [[len(y) for y in x] for x in split]
        max_positions = [x.index(max(x)) for x in lengths]
        res = [x[max_positions[i]] for i, x in enumerate(split)]
        self.needle = [needle(ref, x) for x in res]
        self.alignments = [x[0] for x in self.needle]
        self.scores = [x[1] for x in self.needle]
        #print self.scores
        for i, score in enumerate(self.scores):
            if score < 1300:
                new_needle = needle(ref, reverse_complement(res[i]))
                self.alignments[i] = new_needle[0]
                self.scores[i] = new_needle[1]
        self.needle = []
        for i, alignment in enumerate(self.alignments):
            self.needle.append((alignment, self.scores[i]))

        # Extract aligned sequences
        # Actually, need n ref sequences in case of insertions
        # - ref would have - inserted in that case.
        a_refs = [x[0].seq.tostring().upper() for x in self.alignments]
        a_res = [x[1].seq.tostring().upper() for x in self.alignments]

        for i, reference in enumerate(a_refs):
            frontcount = 0
            backcount = 0
            while True:
                if reference[-backcount - 1] == '-':
                    backcount += 1
                else:
                    break
            while True:
                if reference[frontcount] == '-':
                    frontcount += 1
                else:
                    break
            a_refs[i] = a_refs[i][:len(a_refs[i]) - backcount]
            a_res[i] = a_res[i][:len(a_res[i]) - backcount]
            a_refs[i] = a_refs[i][frontcount:]
            a_res[i] = a_res[i][frontcount:]

        self.ref = a_refs
        self.res = a_res

        self.aligned = [(a_refs[i], a_res[i]) for i, x in enumerate(a_refs)]

        # degapping
        self.gaps = []
        for i, result in enumerate(a_res):
            new_gap = ''.join([_findgap([z]) for j, z in enumerate(result)])
            self.gaps.append(new_gap)
        self.ranges = []
        for gap in self.gaps:
            first = gap.find('X')
            last = (len(gap) - gap[::-1].find('X'))
            self.ranges.append((first, last))

        #mismatches
        self.mismatches = {}
        self.insertions = {}
        self.deletions = {}
        for i, reference in enumerate(a_refs):
            curkey = str(i)
            #curkey = str(i) + ':' + self.resnames[i]
            for j in range(len(reference)):
                if a_refs[i][j] != '-' and a_res[i][j] != '-':
                    if a_refs[i][j] != a_res[i][j]:
                        if curkey in self.mismatches:
                            self.mismatches[curkey].append((j, j))
                        else:
                            self.mismatches[curkey] = [(j, j)]
                elif a_refs[i][j] == '-' and a_res[i][j] != '-':
                    if j in range(self.ranges[i][0], self.ranges[i][1]):
                        if curkey in self.insertions:
                            self.insertions[curkey].append((j, j))
                        else:
                            self.insertions[curkey] = [(j, j)]
                elif a_refs[i][j] != '-' and a_res[i][j] == '-':
                    if j in range(self.ranges[i][0], self.ranges[i][1]):
                        if curkey in self.deletions:
                            self.deletions[curkey].append((j, j))
                        else:
                            self.deletions[curkey] = [(j, j)]

        # Group deletions that happen all in a row.
        for key, value in self.deletions.iteritems():
            last = (-2, -2)
            newlist = []
            for j, deletion in enumerate(value):
                if deletion[0] == last[1] + 1:
                    newlist[-1] = (last[0], deletion[1])
                else:
                    newlist.append(deletion)
                last = newlist[-1]
            self.deletions[key] = newlist

        # Group insertions that happen all in a row.
        for key, value in self.insertions.iteritems():
            last = (-2, -2)
            newlist = []
            for j, deletion in enumerate(value):
                if deletion[0] == last[1] + 1:
                    newlist[-1] = (last[0], deletion[1])
                else:
                    newlist.append(deletion)
                last = newlist[-1]
            self.insertions[key] = newlist

    def report(self):
        '''Report deletions, mismatches, and insertions.'''

        print 'mismatches: ' + str(self.mismatches)
        print 'insertions: ' + str(self.insertions)
        print 'deletions: ' + str(self.deletions)
        if len(self.mismatches) > 0:
            print '----------------------'
            print '----- MISMATCHES -----'
            print '----------------------'
            for key, value in self.mismatches.iteritems():
                index = int(key)
                print '  ' + self.resnames[index]
                for mismatch in value:
                    print ''
                    refi = self.ref[index]
                    resi = self.res[index]
                    _sequences_display([refi, resi], mismatch)
                    print ''
        if len(self.insertions) > 0:
            print '----------------------'
            print '----- INSERTIONS -----'
            print '----------------------'
            for key, value in self.insertions.iteritems():
                index = int(key)
                print '  ' + self.resnames[index]
                for insertion in value:
                    print ''
                    refi = self.ref[index]
                    resi = self.res[index]
                    _sequences_display([refi, resi], insertion)
                    print ''
        if len(self.deletions) > 0:
            print '---------------------'
            print '----- DELETIONS -----'
            print '---------------------'
            for key, value in self.deletions.iteritems():
                index = int(key)
                print '  ' + self.resnames[index]
                for deletion in value:
                    print ''
                    refi = self.ref[index]
                    resi = self.res[index]
                    _sequences_display([refi, resi], deletion)
                    print ''

    def plot(self):
        '''Plot visualization of results using matplotlib.'''

        # Step 1: turn results into ranges, bin those ranges for displaying
        bins = _disjoint_bins(self.ranges)

        # Step 2: plot 'reference' bar
        fig = pyplot.figure()
        sub1 = fig.add_subplot(111)
        gr0 = (0, len(self.aligned[0][0]))
        sub1.broken_barh([gr0], (14.25, 1), facecolors='black',
                         edgecolors='none')

        # Plot and color features
        max_len = len(self.ref_raw.features)

        # Bin the features for plotting
        features = self.ref_raw.features
        # + 1 to start since it seems to be indexed starting with 0
        feature_ranges = [(feature.location.start.position + 1,
                           feature.location.end.position)
                          for feature in features]
        feature_bins = _disjoint_bins(feature_ranges)
        feature_nbin = max(feature_bins)

        for i, feature in enumerate(self.ref_raw.features):
            feature_bin = feature_bins[i]
            qual = feature.qualifiers['label'][0]
            feature_start = feature.location.start.position
            feature_end = feature.location.end.position
            mid = int(math.ceil((feature_start + feature_end) / 2))
            centered = (feature_bin + 1) * 10
            sub1.broken_barh([(feature_start, feature_end-feature_start)],
                             (centered, 9),
                             facecolors=cm.Set3(float(i) / max_len),
                             edgecolors='black')
            sub1.text(mid, centered + 7, qual,
                      rotation=90)

        def add_discrepancies(index, height):
            '''
            Add insertions, deletions, and mismatches to plot.

            :param height: height of the annotation (on plot).
            :type bin_index: int

            '''

            sub1.plot(1000, 25)
            index = str(index)

            for key, value in self.insertions.iteritems():
                if key == index:
                    for insertion in value:
                        sub1.plot(insertion[0], height, marker='o', color='k')
            for key, value in self.deletions.iteritems():
                if key == index:
                    for deletion in value:
                        sub1.plot(deletion[0], height, marker='^', color='k')
            for key, value in self.mismatches.iteritems():
                if key == index:
                    for mismatch in value:
                        sub1.plot(mismatch[0], height, marker='*', color='k')

        def wrap_name(str_in, wrap_len=10):
            '''
            Wrap plotted text to avoid overlaps (not perfect).

            :param str_in: Input string.
            :type str_in: str
            :param wrap_len: Length at which to wrap lines.
            :type wrap_len: str

            '''

            wrap_positions = range(0, len(str_in), wrap_len)
            out = [str_in[i:i + wrap_len] for i in wrap_positions]
            return '\n'.join(out)

        # Step 3: plot results ranges
        for i, current_range in enumerate(self.ranges):
            for vals in current_range:
                gap = self.gaps[i]
                gap_index = gap.find('X')
                ends = (gap_index, len(gap) - gap[::-1].find('X') - gap_index)
                centered = (bins[i] + 1) * 10 + 10 + 10 * feature_nbin
                sub1.broken_barh([ends], (centered, 9), facecolors='pink',
                                 edgecolors='black')
                sub1.text(ends[0] + 10, centered + 8,
                          wrap_name(self.resnames[i]),
                          verticalalignment='top')
                add_discrepancies(i, centered + 2)

        sub1.set_xlim(0, gr0[1])
        sub1.set_xlabel('Base pairs from origin')
        sub1.set_yticks([15, 15 + 10 * (feature_nbin + 1)])
        sub1.set_yticklabels(['Reference', 'Results'])
        sub1.grid(True)
        sub1.xaxis.grid(False)
        sub1.yaxis.grid(False)

        pyplot.title('Alignment gap summary', fontstyle='normal')
        pyplot.show()

    def write_alignment(self):
        '''
        Write alignment results to file - allows reanalysis and processing
        by other programs.

        '''

        # custom format or standard (e.g. FASTA)? Implement both?
        pass

    def write_plot(self):
        '''Save plot to image (png or svg).'''

        pass

    def fix_disrepancy(self, result, position, newvalue):
        '''Fix mismatch, deletion, or insertion manually.'''

        pass


def readref(filepath, ftype='genbank'):
    '''
    Read reference genbank file.

    :param filepath: Path to the file (typically genbank).
    :type filepath: str
    :param ftype: Valid filetype readable by Bio.SeqIO.
    :type ftype: str

    '''

    sequence = SeqIO.read(filepath, ftype)
    return sequence


def readres(dirpath):
    '''
    Read .seq results files from a dir.

    :param dirpath: Path to directory containing sequencing files.
    :type dirtpath: str

    '''

    seq_paths = [x for x in os.listdir(dirpath) if x.endswith('.seq')]
    abi_paths = [x for x in os.listdir(dirpath) if x.endswith('.abi')]
    abi_paths += [x for x in os.listdir(dirpath) if x.endswith('.ab1')]
    seq_seqs = [SeqIO.read(dirpath + x, 'fasta') for x in seq_paths]
    abi_seqs = [SeqIO.read(dirpath + x, 'abi') for x in abi_paths]
    sequences = seq_seqs + abi_seqs
    return sequences


def _findgap(list_in):
    '''
    Iterate over string list, return 'X' if has non-\'-\' value.

    :param list_in: String list.
    :type list_in: list

    '''

    for base in list_in:
        if base != '-':
            return 'X'
    return '-'


def _sequences_display(seqs, start_stop, context=10, indent=4):
    '''
    Given two sequences to compare, display them and visualize non-matching
    regions.

    :param seqs: Sequences to display.
    :type seqs: list
    :param start_stop: Indices to display.
    :type start_stop: tuple
    :param context: Extra context to add on either side of the displayed
                    sequences.
    :type context: int
    :param indent: Indentation of the displayed text.
    :type indent: int

    '''

    if len(seqs) != 2:
        raise ValueError('Expected two sequences')
    start = start_stop[0]
    stop = start_stop[1]

    seq_list_1 = []
    seq_list_2 = []
    seq_list_1.append(seqs[0][max(start - context, 0):start])
    seq_list_2.append(seqs[1][max(start - context, 0):start])
    seq_list_1.append(seqs[0][start:stop + 1])
    seq_list_2.append(seqs[1][start:stop + 1])
    seq_list_1.append(seqs[0][stop + 1:min(stop + context, len(seqs[0]))])
    seq_list_2.append(seqs[1][stop + 1:min(stop + context, len(seqs[0]))])

    # No bars for non-matching sequences
    top = ''.join(seq_list_1)
    bottom = ''.join(seq_list_2)
    middle = ['|' if top[i] == bottom[i] else ' ' for i in range(len(top))]
    middle = ''.join(middle)

    indent = ' ' * indent
    if start != stop:
        print indent + 'Positions {0} to {1}:'.format(start, stop)
    else:
        print indent + 'Position {}:'.format(start)
    print indent + top
    print indent + middle
    print indent + bottom


def _disjoint_bins(range_tuple_list):
    '''
    Construct disjoint bins given a list of range tuples.

    :param range_tuple_list: A list of tuples containing range values.
    :type range_tuple_list: list

    '''

    rtl = range_tuple_list
    # number the ranges (third value in tuple)
    rtl = [(x[0], x[1], i) for i, x in enumerate(rtl)]

    # sort by range start
    rtl = sorted(rtl, key=lambda starts: starts[0])
    rtl_len = len(rtl)
    remaining = rtl

    done_binning = False
    binned = []
    while not done_binning:
        current_bin = []
        nextbin = []
        while len(remaining) > 0:
            last = remaining.pop(0)
            current_bin.append(last)
            nextbin += [x for x in remaining if x[0] <= last[1]]
            remaining = [x for x in remaining if x[0] > last[1]]
        binned.append(current_bin)
        if len(remaining) == 0 and len(nextbin) == 0:
            done_binning = True
        else:
            remaining = nextbin

    bin_list = [0] * rtl_len
    for i, bin_ranges in enumerate(binned):
        for bin_range in bin_ranges:
            bin_list[bin_range[2]] = i

    return bin_list
