'''
Sanger sequencing alignment tools.

'''

import math
from matplotlib import pyplot
from matplotlib import cm

from pymbt.analysis._sequencing.needle import needle

# TODO:
# consensus / master sequence for plotting / report / analysis
# too much work is done initializing Sanger. Move calculations to methods.


class Sanger(object):
    '''
    Align and analyze Sanger sequencing results.

    '''

    def __init__(self, reference, results):
        '''
        :param reference: Reference sequence.
        :type reference: DNA object
        :param results: Sequencing result string. A list of DNA objects is also
                        valid.
        :type results: DNA object

        '''

        # Make results a list for consistency
        if type(results) != list:
            results = [results]

        # Coerce inputs to string
        results = [str(seq) for seq in results]
        for seq in results:
            seq = str(seq)
        reference = str(reference)

        self.ref_raw = reference
        self.resnames = [x.name for x in results]

        # Reduce results to largest non-N segment
        for seq in results:
            split = seq.split('N')
            sizes = [len(x) for x in split]
            seq = split[sizes.index(max(sizes))]

        # Align
        self.needle = [needle(reference, result) for result in results]
        self.alignments = self.needle['alignments']
        self.scores = self.needle['scores']

        # If score is too low, try reverse complement of result
        for i, score in enumerate(self.scores):
            if score < 1300:
                new_needle = needle(reference, results[i].reverse_complement())
                self.alignments[i] = new_needle['alignments']
                score = new_needle['scores']

        # Extract aligned sequences
        # Actually, need n ref sequences in case of insertions
        # - ref would have - inserted in that case.
        aligned_refs = [x[0].seq.tostring().upper() for x in self.alignments]
        aligned_res = [x[1].seq.tostring().upper() for x in self.alignments]

        for a_ref, a_reses in zip(aligned_refs, aligned_res):
            frontcount = 0
            backcount = 0
            while True:
                if a_ref[-backcount - 1] == '-':
                    backcount += 1
                else:
                    break
            while True:
                if a_ref[frontcount] == '-':
                    frontcount += 1
                else:
                    break
            a_ref = a_ref[:len(a_ref) - backcount]
            a_reses = a_reses[:len(a_reses) - backcount]
            a_ref = a_ref[frontcount:]
            a_reses = a_reses[frontcount:]

        self.ref = aligned_refs
        self.res = aligned_res

        self.aligned = [(a_ref, a_reses) for a_ref, a_reses in
                        zip(aligned_refs, aligned_res)]

        # Degap
        self.gaps = []
        for result in enumerate(aligned_res):
            new_gap = ''.join([_findgap([x]) for x in result])
            self.gaps.append(new_gap)
        self.ranges = []
        for gap in self.gaps:
            first = gap.find('X')
            last = (len(gap) - gap[::-1].find('X'))
            self.ranges.append((first, last))

        # Mismatches, Insertions, and Deletions
        self.mismatches = {}
        self.insertions = {}
        self.deletions = {}
        # TODO: This is unreadable. Rewrite it.
        for i, reference in enumerate(aligned_refs):
            i_key = str(i)
            for j in range(len(reference)):
                if aligned_refs[i][j] != '-' and aligned_res[i][j] != '-':
                    if aligned_refs[i][j] != aligned_res[i][j]:
                        if i_key in self.mismatches:
                            self.mismatches[i_key].append((j, j))
                        else:
                            self.mismatches[i_key] = [(j, j)]
                elif aligned_refs[i][j] == '-' and aligned_res[i][j] != '-':
                    if j in range(self.ranges[i][0], self.ranges[i][1]):
                        if i_key in self.insertions:
                            self.insertions[i_key].append((j, j))
                        else:
                            self.insertions[i_key] = [(j, j)]
                elif aligned_refs[i][j] != '-' and aligned_res[i][j] == '-':
                    if j in range(self.ranges[i][0], self.ranges[i][1]):
                        if i_key in self.deletions:
                            self.deletions[i_key].append((j, j))
                        else:
                            self.deletions[i_key] = [(j, j)]

        # Group deletions that happen all in a row.
        for key, value in self.deletions.iteritems():
            last = (-2, -2)
            newlist = []
            for deletion in value:
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
        '''
        Report deletions, mismatches, and insertions.

        '''

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
                    _sequences_display([refi, resi], mismatch[0], mismatch[1])
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
                    _sequences_display([refi, resi], insertion[0],
                                       insertion[1])
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
                    _sequences_display([refi, resi], deletion[0], deletion[1])
                    print ''

    def plot(self):
        '''
        Plot visualization of results using matplotlib.

        '''

        # FIXME: have to implement features in order to plot them

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
            :type height: int

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
        '''
        Save plot to image (png or svg).

        '''

        pass

    def fix_disrepancy(self, result, position, newvalue):
        '''
        Fix mismatch, deletion, or insertion manually.

        '''

        pass


# TODO: why is this necessary?
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


def _sequences_display(seqs, start, stop, context=10, indent=4):
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

    seq_lists = [[], []]

    for seq, seq_list in zip(seqs, seq_lists):
        seq_list.append(seq[max(start - context, 0):start])
        seq_list.append(seq[start:stop + 1])
        seq_list.append(seq[stop + 1:min(stop + context, len(seq))])

    # No bars for non-matching sequences
    top = ''.join(seq_list[0])
    bottom = ''.join(seq_list[1])
    middle = ['|' if t == b else ' ' for t, b in zip(top, bottom)]
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

    remaining = rtl[:]

    binned = []
    while True:
        current_bin = []
        next_bin = []
        while remaining:
            last = remaining.pop(0)
            current_bin.append(last)
            next_bin += [x for x in remaining if x[0] <= last[1]]
            remaining = [x for x in remaining if x[0] > last[1]]
        binned.append(current_bin)

        if not remaining and not next_bin:
            break
        else:
            remaining = next_bin

    bin_list = [0] * rtl_len
    for i, bin_ranges in enumerate(binned):
        for bin_range in bin_ranges:
            bin_list[bin_range[2]] = i

    return bin_list
