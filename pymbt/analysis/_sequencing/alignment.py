'''
Sanger sequencing alignment tools.

'''

from matplotlib import pyplot
from matplotlib import cm

from pymbt.analysis import needle

# TODO:
# consensus / master sequence for plotting / report / analysis


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

        self.reference = reference
        self.results_input = results
        self.names = [seq.name for seq in results]

        # Remove Ns from sequencing results - keep largest remaining segment
        self.processed = [self._remove_n(x) for x in results]

        # Align
        print '(Aligning...)'
        self.alignments, self.scores = self._align()

        # Find leading/trailing '-' and record coverage indices
        self.coverage = self._find_coverage()

        # Calculate mismatches, insertions, and deletions
        discrepancies = self._find_discrepancies()
        self.mismatches = discrepancies[0]
        self.insertions = discrepancies[1]
        self.deletions = discrepancies[2]

    def report(self):
        '''
        Report deletions, mismatches, and insertions.

        '''

        n_mismatches = sum([len(x) for x in self.mismatches])
        n_insertions = sum([len(x) for x in self.insertions])
        n_deletions = sum([len(x) for x in self.deletions])

        mismatch_dict = {'number': n_mismatches,
                         'name': 'mismatches',
                         'list': self.mismatches}
        insertions_dict = {'number': n_insertions,
                           'name': 'insertions',
                           'list': self.insertions}
        deletions_dict = {'number': n_deletions,
                          'name': 'deletions',
                          'list': self.deletions}
        disparities = [mismatch_dict, insertions_dict, deletions_dict]

        col1 = 'mismatches: {}  '.format(n_mismatches)
        col2 = 'insertions: {}  '.format(n_insertions)
        col3 = 'deletions: {}'.format(n_deletions)
        print
        print 'Report: '
        print
        print col1 + col2 + col3

        for disparity in disparities:
            if disparity['number']:
                print '----------------------'
                print '----- {}'.format(disparity['name'])
                print '----------------------'

                for i, result_disparity in enumerate(disparity['list']):
                    result_name = self.names[i]
                    print '  {}'.format(result_name)
                    ref_i = self.alignments[i][0]
                    res_i = self.alignments[i][1]

                    for start, end in result_disparity:
                        print
                        _sequences_display(ref_i, res_i, start, end)
                        print

    def plot(self):
        '''
        Plot visualization of the alignment results using matplotlib.

        '''

        # Plots:
        # 1. reference bar at the bottom
        # 2. binned, labeled reference sequence features on top of that
        # 3. binned, labeled sequencing results above reference features
        # 4. discrepancies on top of sequencing results

        ################
        # Calculations #
        ################

        # Reference bar information
        reference_x = 0
        reference_y = 14.25
        reference_width = len(self.alignments[0][0])
        reference_height = 1

        # Bin the features so they don't overlap when plotted
        features = self.reference.features
        feature_ranges = [(feature.start, feature.stop) for feature in
                          features]
        feature_bins = _disjoint_bins(feature_ranges)
        feature_nbin = len(feature_bins)

        # Bin the alignments so they don't overlap when plotted
        alignment_bins = _disjoint_bins(self.coverage)

        # Calculate discrepancy coordinates
        discrepancy_coords = [[], [], []]
        discrepancies = [self.mismatches, self.insertions, self.deletions]
        for i, result_bin in enumerate(alignment_bins):
            for index in result_bin:
                for j, discrepancy_type in enumerate(discrepancies):
                    for discrepancy in discrepancy_type[index]:
                        y = i
                        x = (discrepancy[0] + discrepancy[1]) // 2
                        discrepancy_coords[j].append((x, y))

        ############
        # Plotting #
        ############

        # Plot calculations
        # Controls spacing between bars, size of bars, etc
        size = 10

        # Matplotlib commands
        # Plot a black reference bar at the bottom
        fig = pyplot.figure()
        sub1 = fig.add_subplot(111)
        sub1.broken_barh([(reference_x, reference_width)],
                         (reference_y, reference_height),
                         facecolors='black', edgecolors='none')

        # Plot the reference features on top of the bar
        # Plot the features by bin:
        for i, feature_bin in enumerate(feature_bins):
            for index in feature_bin:
                feature = features[index]
                name = feature.name

                width = feature.stop - feature.start
                height = size - 1
                y_index = (i + 1) * size
                mid = (feature.start + feature.stop) // 2

                pos = float(i) / len(features)
                sub1.broken_barh([(feature.start, width)],
                                 (y_index, height),
                                 facecolors=cm.Set3(pos),
                                 edgecolors='black')
                sub1.text(mid, y_index + size / 2, name, rotation=90)

        # Plot sequencing results by bin
        for i, result_bin in enumerate(alignment_bins):
            for index in result_bin:
                start, stop = self.coverage[index]
                name = self.names[index]

                width = stop - start
                height = size - 1
                y_index = (i + 1) * size + size * (feature_nbin + 1)
                text_x = start + (stop - start) // 6

                sub1.broken_barh([(start, width)], (y_index, height),
                                 facecolors='pink', edgecolors='black')
                sub1.text(text_x, y_index + size // 3, _wrap_name(name),
                          rotation=0)

        # Plot mismatches, insertions, deletions
        sub1.plot(1000, 25)
        shapes = ['o', 'x', 'v']
        labels = ['mismatch', 'insertion', 'deletion']
        for coords, shape, label in zip(discrepancy_coords, shapes, labels):
            x_coords = [x[0] for x in coords]
            y_coords = [(x[1] + feature_nbin + 2) * size for x in coords]
            sub1.scatter(x_coords, y_coords, marker=shape, color='k',
                         label=label)
        sub1.legend()

        # Plot labeling, etc
        sub1.set_xlim(0, reference_width)
        sub1.set_xlabel('Base pairs from origin')
        sub1.set_yticks([15, 15 + size * (feature_nbin + 1)])
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
        return NotImplemented

    def write_plot(self):
        '''
        Save plot to image (png or svg).

        '''

        return NotImplemented

    def fix_disrepancy(self, result, position, newvalue):
        '''
        Fix mismatch, deletion, or insertion manually.

        '''

        return NotImplemented

    def _remove_n(self, seq):
        '''
        Remove Ns from sequence - strategy is to find largest non-N segment
        and return it.

        :param seq: Sequence that contains Ns to remove
        :type seq: str
        '''

        split = seq.top.split('n')
        sizes = [len(x) for x in split]
        largest = split[sizes.index(max(sizes))]
        seq_start = seq.top.index(largest)
        seq_stop = seq_start + len(largest)
        processed = seq[seq_start:seq_stop]

        return processed

    def _align(self):
        '''
        Aligns sequences in Sanger. self.reference and self.processed have to
        exist first.

        '''

        # Align
        needle_result = [needle(str(self.reference), str(seq)) for seq in
                         self.processed]

        # Split into alignments and scores
        alignments = [(result[0], result[1]) for result in needle_result]
        scores = [result[2] for result in needle_result]

        # If score is too low, may be a 'reverse' sequencing reaction
        # Try reverse complement
        for i, score in enumerate(scores):
            if score < 1300:
                reversed_result = self.processed[i].reverse_complement()
                new_needle = needle(str(self.reference),
                                    str(reversed_result))
                alignments[i] = (new_needle[0], new_needle[1])
                score = new_needle[2]

        return alignments, scores

    def _find_coverage(self):
        '''
        Finds coverage of the alignments, i.e. figures out where trailing
        and leading gaps are. self.alignments has to exist.

        '''

        coverage = []
        for ref, res in self.alignments:
            frontcount = 0
            backcount = 0
            while res[-backcount - 1] == '-':
                backcount += 1
            while res[frontcount] == '-':
                frontcount += 1
            coverage.append((frontcount, len(ref) - backcount))

        return coverage

    def _find_discrepancies(self):
        '''
        Find discrepancies in the alignments - mismatches, insertions, and
        deletions.

        self.alignments and self.coverage must exist for this to work.

        '''

        mismatches = [[] for i in range(len(self.alignments))]
        insertions = [[] for i in range(len(self.alignments))]
        deletions = [[] for i in range(len(self.alignments))]
        # Got through each alignment
        for i, (reference, result) in enumerate(self.alignments):
            # Go through each alignment base by base
            coverage = self.coverage[i]
            ref_trim = reference[slice(*coverage)]
            res_trim = result[slice(*coverage)]
            for j, (refbase, resbase) in enumerate(zip(ref_trim, res_trim)):
                # If only the reference base is a '-', is insertion
                position = j + coverage[0]
                if refbase == '-':
                    insertions[i].append((position, position))
                # If only the reference base is a '-', is deletion
                elif resbase == '-':
                    deletions[i].append((position, position))
                # If bases don't match, is mismatch
                elif refbase != resbase:
                    mismatches[i].append((position, position))

        # Group mismatches/indels if they appear one after another
        for i, alignment_mismatches in enumerate(mismatches):
            mismatches[i] = _group_indels(alignment_mismatches)
        for i, alignment_insertions in enumerate(insertions):
            insertions[i] = _group_indels(alignment_insertions)
        for i, alignment_deletions in enumerate(deletions):
            deletions[i] = _group_indels(alignment_deletions)

        return mismatches, insertions, deletions


def _group_indels(indel_list):
    '''
    Group insertions or deletions that form a continuous block.

    :param indel_list: list of deletion or insert indices (elements are
                       2-tuples)
    :type indel_list: list
    '''

    # By setting to (-2, -2), skips first base so that others can
    # 'look back' to see whether the current deletion is last + 1
    previous = (-2, -2)
    grouped_indels = []
    for indel in indel_list:
        if indel[0] - previous[1] == 1:
            grouped_indels[-1] = (previous[0], indel[1])
        else:
            grouped_indels.append(indel)
        previous = grouped_indels[-1]

    return grouped_indels


def _sequences_display(seq1, seq2, start, stop, context=10, indent=4):
    '''
    Given two sequences to compare, display them and visualize non-matching
    regions. Sequences should be about the same size, ideally.

    :param seq1: First sequence to compare.
    :type seq1: str
    :param seq2: Second sequence to compare.
    :type seq2: str
    :param start: Where to start displaying the sequence.
    :type start: int
    :param stop: Where to stop displaying the sequence.
    :type stop: int
    :param context: Extra context to add on either side of the displayed
                    sequences.
    :type context: int
    :param indent: Indentation of the displayed text.
    :type indent: int

    '''
    # TODO: if seqs aren't the same size, should display a little differently
    # TODO: display is incomplete - doesn't show all mismatches / deletions /
    # insertions at once, just the ones displayed at the moment.
    # Should instead calculate a single 'overview' set of sequences:
    # top, bottom, and middle, where middle is where all the seqs match.
    # Then just subset this

    # Figure out how much context can be included (seq might start/end earlier)
    l_context = (max(start - context - 1, 0), start - 1)
    r_context = (stop + 1, min(start + context + 1, len(seq1)))

    def gen_levels(seq1, seq2, context_tuple):
        '''
        Generate a column to display - left side or right side of disparity.

        :param seq1: same as seq1 of parent
        :type seq1: str
        :param seq2: same as seq2 of parent
        :type seq2: str
        :context_tuple: l_context or r_context
        :type context_tuple: tuple

        '''

        top, bottom = [seq[slice(*context_tuple)] for seq in [seq1, seq2]]
        context_len = (context_tuple[1] - context_tuple[0])
        middle = '|' * context_len
        highlight = ' ' * context_len
        return top, bottom, middle, highlight

    # Generate left and right columns of text to display
    l_top, l_bottom, l_middle, l_highlight = gen_levels(seq1, seq2, l_context)
    r_top, r_bottom, r_middle, r_highlight = gen_levels(seq1, seq2, r_context)

    # Generate core column - core is where discrepancies are
    core_slice = slice(start, stop + 1)
    core_len = stop - start + 1
    core_top = seq1[core_slice]
    core_bottom = seq2[core_slice]
    core_middle = ' ' * core_len
    core_highlight = '*' * core_len

    # Combine each level together for printing
    top = l_top + core_top + r_top
    middle = l_middle + core_middle + r_middle
    bottom = l_bottom + core_bottom + r_bottom
    highlight = l_highlight + core_highlight + r_highlight

    indent = ' ' * indent
    if start != stop:
        print '{0}Positions {1} to {2}:'.format(indent, start, stop)
    else:
        print '{0}Position {1}:'.format(indent, start)

    print indent + top
    print indent + middle
    print indent + bottom
    print indent + highlight


def _disjoint_bins(ranges_list):
    '''
    Construct disjoint bins given a list of 1-D ranges (tuples). Returns a list
    of bins containing the index of each range that fell into that bin.

    :param range_tuple_list: A list of tuples containing range values.
    :type range_tuple_list: list

    '''

    # Keep track of the original order for reporting later
    ranges_list = [(x[0], x[1], i) for i, x in enumerate(ranges_list)]
    ranges_list = sorted(ranges_list, key=lambda starts: starts[0])

    remaining = ranges_list[:]
    binned = []
    while True:
        current_bin = []
        next_bin = []
        while remaining:
            last = remaining.pop(0)
            current_bin.append(last)
            next_bin += [x for x in remaining if x[0] < last[1]]
            remaining = [x for x in remaining if x[0] >= last[1]]
        binned.append(current_bin)

        if not remaining and not next_bin:
            break
        else:
            remaining = next_bin

    bins_pre = [[x[2] for x in range_bin] for range_bin in binned]
    bins = bins_pre

    return bins


def _wrap_name(str_in, wrap_len=15):
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
