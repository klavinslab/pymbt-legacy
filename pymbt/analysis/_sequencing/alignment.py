'''
Sanger sequencing alignment tools.

'''

from matplotlib import pyplot
from matplotlib import cm

from pymbt.analysis._sequencing.needle import needle

# bigger TODO: This module doesn't work at all right now. Fix it.
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
        :type reference: pymbt.sequence.DNA
        :param results: Sequencing result string. A list of DNA objects is also
                        valid.
        :type results: list of pymbt.sequence.DNA elements

        '''

        # Make results a list for consistency
        if type(results) != list:
            results = [results]

        self.reference_input = reference
        self.results_input = results

        # Remove Ns from sequencing results - keep largest remaining segment
        self.processed = []
        for seq in results:
            split = seq.top.split('n')
            sizes = [len(x) for x in split]
            largest = split[sizes.index(max(sizes))]
            seq_start = seq.top.index(largest)
            seq_stop = seq_start + len(largest)
            self.processed.append(seq[seq_start:seq_stop])

        # Align
        needle_result = [needle(str(reference), str(seq)) for seq in
                         self.processed]
        self.alignments = [result['alignments'][0] for result in needle_result]
        self.scores = [result['scores'][0] for result in needle_result]

        # If score is too low, may be a 'reverse' sequencing reaction
        # Try reverse complement
        for i, score in enumerate(self.scores):
            if score < 1300:
                new_needle = needle(str(reference),
                                    str(results[i].reverse_complement()))
                self.alignments[i] = new_needle['alignments'][0]
                score = new_needle['scores']

        # Find leading/trailing '-' and 1) note their location and 2)
        # make versions that remove them.
        self.alignment_coverage = []
        for ref, res in zip(*self.alignments):
            frontcount = 0
            backcount = 0
            while ref[-backcount - 1] == '-':
                backcount += 1
            while ref[frontcount] == '-':
                frontcount += 1
            frontcount
            len(ref) - backcount
            len(res) - backcount

        # Find gaps in sequencing results
        self.gaps = []
        for reference, result in zip(*self.alignments):
            temp_gap = []
            for i, ref_base in enumerate(reference):
                if ref_base == '-' and results[i] == '-':
                    temp_gap.append(i)
            self.gaps.append(temp_gap)


        # Based on gaps, find ranges for which there is coverage (no gaps)
        self.ranges = []
        for gap in self.gaps:
            first = gap.find('X')
            last = (len(gap) - gap[::-1].find('X'))
            self.ranges.append((first, last))

        # Calculate mismatches, insertions, and deletions
        self.mismatches = []
        self.insertions = []
        self.deletions = []

        # TODO: This is unreadable. Rewrite it.
        # Got through each alignment
        for ref, res in zip(*self.alignments):
            # Go through each alignment base by base
            for j, refbase in enumerate(ref):
                # If both bases are a gap, it's a terminal gap (ignore it)
                if refbase == '-' and res[j] == '-':
                    pass
                # If only the reference base is a '-', report insertion
                elif refbase == '-':
                    if j in range(self.ranges[i][0], self.ranges[i][1]):
                        if i_key in self.insertions:
                            self.insertions[i_key].append((j, j))
                        else:
                            self.insertions[i_key] = [(j, j)]
                # If only the reference base is a '-', report deletion
                elif res[j] == '-':
                    if j in range(self.ranges[i][0], self.ranges[i][1]):
                        if i_key in self.deletions:
                            self.deletions[i_key].append((j, j))
                        else:
                            self.deletions[i_key] = [(j, j)]
                # If bases don't match, report mismatch
                elif refbase != res[j]:
                    if i_key in self.mismatches:
                        self.mismatches[i_key].append((j, j))
                    else:
                        self.mismatches[i_key] = [(j, j)]

        # Deletions and insertions are found on a per-base basis. If they occur
        # one after another, group together as a single deletion or insertion
        # Deletion grouping:
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
        # Insertion grouping:
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

        print 'mismatches: {}'.format(self.mismatches)
        print 'insertions: {}'.format(self.insertions)
        print 'deletions: {}'.format(self.deletions)
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
                print '  {}'.format(self.resnames[index])
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
                print '  {}'.format(self.resnames[index])
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
            mid = (feature_start + feature_end) // 2
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


def _findgap(result):
    '''
    Iterate over string list, return 'X' if has non-\'-\' value.

    :param result: String list.
    :type result: list

    '''

    for base in result:
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
        print '{0}Positions {1} to {2}:'.format(indent, start, stop)
    else:
        print '{0}Position {1}:'.format(indent, start)
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
