'''Split genes or arbitrary DNA sequences into chunks for optimal assembly.'''

import csv
import itertools
import time
from pymbt.analysis import StructureWindows

# TODO: circular (e.g. plasmid) version.
#   1. Convert plasmid to 'linear' version - add 60bp from front to end
#   2. Evaluate as if linear
#   3. Rotate sequence by the window size and repeat until whole plasmid
#   has been analyzed
#   4. Pick the best one
#   Alternative:
#       Since we'd be analyzing the whole sequences anyways, just do that
#       right off the bat, store it, and try making spanning combinations


class GeneSplitter(object):
    '''
    Split a large (~<10kb) sequence into smaller ones that are easier to clone,
    either via oligo assembly or PCR. Store and write the results.

    '''

    def __init__(self, dna_object):
        '''
        :param seq: Input sequence - a DNA object.
        :type seq: DNA object

        '''

        self.template = dna_object
        self.overlaps = []
        self.scores = []
        self.split_sequences = []

    def split(self, max_len=1100, min_context=200, core=60, context=90,
              step=10, force_exhaustive=False):
        '''
        Split a large (~<10kb) sequence into smaller ones that are easier to
        clone, either via oligo assembly or PCR.

        :param max_len: Maximum length of each split-up output sequence.
        :type max_len: int
        :param min_context: Sets minimum size of each potential overlap
                            region to evaluate.
        :type min_context: int
        :param core: Window size.
        :type core: int
        :param context: Amount (in bp) of context to evaluate on either side
                        of the window.
        :type context: int
        :param step: Step size for the windows.
        :type step: int
        :param force_exhaustive: Forces algorithm to find global optimum
                                 exhaustively. Can dramatically slow down the
                                 computation.
        :type force_exhaustive: bool

        '''

        seq_len = len(self.template)
        overlap_n = 1
        while True:
            if seq_len < overlap_n * max_len - (overlap_n - 1) * core:
                break
            else:
                overlap_n += 1

        # calculate potential overlap regions
        if overlap_n == 1:
            # doesn't need fancy calculations - no overlaps
            return self.template
        elif overlap_n > 1:
            while True:
                # Generate overlap indices.
                # If overlaps smaller than min_context, increment n
                # Incrementing n in this way can lead to huge areas
                # for potential overlaps - may want to trim down for
                # redundant nupack usage
                stops = []
                for i in range(overlap_n - 1):
                    stops.append((i + 1) * max_len - i * core)
                stops = [x if x < seq_len else seq_len for x in stops]
                starts = [seq_len - x for x in stops]
                starts = [x if x > 0 else 0 for x in starts]
                starts.reverse()
                olap_start_stop = [[starts[i], x] for i, x in enumerate(stops)]
                olaps = [self.template[x[0]:x[1]] for x in olap_start_stop]
                with_context = []
                for ends in olap_start_stop:
                    start_i = max(ends[0] - context, 0)
                    stop_i = min(ends[1] + context, seq_len)
                    with_context.append([start_i, stop_i, seq_len])
                olaps_w_context = [self.template[x[0]:x[1]] for x in
                                   with_context]
                olap_w_context_lens = [x[1] - x[0] for x in with_context]
                if not any([x < min_context for x in olap_w_context_lens]):
                    break
                overlap_n += 1
        else:
            raise ValueError('Number of pieces is somehow zero.\
                             That shouldn\'t happen.')

        walked_raw = []
        for i, overlap in enumerate(olaps_w_context):
            print 'Analyzing {0} of {1} overlap(s).'.format(i + 1, len(olaps))
            walker = StructureWindows(overlap)
            window = walker.run(core_len=core, context_len=context, step=step)
            walked_raw.append(window)

        walked = []
        for i, window_raw in enumerate(walked_raw):
            new_window = []
            for w_vals in window_raw:
                w_vals0 = w_vals[0] + with_context[i][0]
                w_vals1 = w_vals[1] + with_context[i][0]
                w_vals2 = w_vals[2]
                new_window.append((w_vals0, w_vals1, w_vals2))
            walked.append(new_window)

        scores = [[y[2] for y in x] for x in walked]

        best_overlaps = optimal_overlap(walked, max_len,
                                        force_exhaustive=force_exhaustive)
        olap_starts = [x[0] for x in best_overlaps]
        olap_stops = [x[1] for x in best_overlaps]
        scores = [x[2] for x in best_overlaps]
        olap_start_stop = [(x, olap_stops[i]) for i, x in
                           enumerate(olap_starts)]

        best_starts = [0] + [x[0] for x in best_overlaps]
        best_stops = [x[1] for x in best_overlaps] + [seq_len]
        final_seqs = [self.template[best_starts[i]:best_stops[i]] for i in
                      range(len(best_starts))]
        self.split_sequences = final_seqs
        self.overlaps = olap_start_stop
        self.scores = scores

        return {'split_sequences': self.split_sequences,
                'overlaps': self.overlaps,
                'scores': self.scores}

    def write(self, path):
        '''
        Write out results to csv.

        :param path: Full path of the output csv file (including extension).
        :type path: str

        '''

        handle = open(path, 'w')
        writer = csv.writer(handle, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['sequence', 'olap_start', 'olap_stop', 'olap_score'])
        sta = [x[0] for x in self.overlaps] + ['NA']
        sto = [x[1] for x in self.overlaps] + ['NA']
        for i, sequence in enumerate(self.split_sequences):
            writer.writerow([sequence, sta[i], sto[i], self.scores[i]])
        handle.close()


def optimal_overlap(walked, max_len, force_exhaustive=False):
    '''
    Search for the optimal overlap combination.

    :param walked: list of lists of the output to StructureWindows,
                   representing a collection of each potential overlap's
                   evaluated windows.
    :type walked: list
    :param max_len: Maximum length of each split-up output sequence.
    :type max_len: int
    :param force_exhaustive: Forces algorithm to find global optimum
                             exhaustively. Can dramatically slow down the
                             computation.
    :type force_exhaustive: bool

    '''

    exhaustive = True

    # sort results by score
    def sort_tuple(tuple_in):
        '''
        Sort tuple of results by score.

        :param tuple_in: input tuple to be sorted by third value.
        :type tuple_in: tuple

        '''

        return sorted(tuple_in, key=lambda score: score[2], reverse=True)

    def spannable(start_stop_list):
        '''
        Test whether sequence can be spanned given a list of starts and stops
        for the overlaps.

        :param start_stop_list: a list of start and stop position tuples.
        :type start_stop_list: list

        '''

        starts_list = [[y[0] for y in x] for x in start_stop_list]
        stops_list = [[y[1] for y in x] for x in start_stop_list]
        for starts in starts_list:
            starts.sort()
        for stops in stops_list:
            stops.sort()

        # Start at farthest-right start position

        # Get farthest-right start positions
        right_starts = [x[-1] for x in starts_list]
        # Exclude 'start' of last overlap - nothing for it to span to
        right_starts.pop()

        # Get farthest-left stop positions
        left_stops = [x[0] for x in stops_list]
        # Exclude 'stop' of first overlap - nothing for it to span to
        left_stops.pop(0)

        # See whether the difference between any is less than max_len
        diff = [(left_stops[i] - x > max_len) for i, x in
                enumerate(right_starts)]
        if any(diff):
            return False
        return True

#        current_start = starts_list[0][-1]
#        for i, starts in enumerate(starts_list):
#            working_stops = []
#            next_starts = []
#            for j, stop in enumerate(stops_list[i]):
#                if (stop - current_start) <= max_len:
#                    working_stops.append(stop)
#                    next_starts.append(starts[j])
#            if not working_stops:
#                return False
#            else:
#                current_start = next_starts[-1]
#        return True

    def check_useable(tuple_in):
        '''
        Test whether a set of coordinates is smaller than the max_len.

        :param tuple_in: Coordinates tuple - start of first, end of last.
        :type tuple_in: tuple

        '''

        for i in range(len(tuple_in) - 1):
            if tuple_in[i + 1][1] - tuple_in[i][0] > max_len:
                return False
        return True

    sorted_walked = [sort_tuple(x) for i, x in enumerate(walked)]

    position = 1
    while exhaustive:
        current = [x[0:position] for x in sorted_walked]
        if not spannable(current):
            position += 1
        else:
            remove_nonspanning(current, max_len)
            n_combos = 1
            for combo in current:
                n_combos *= len(combo)
            if n_combos > 10000000 and not force_exhaustive:
                exhaustive = False
                break
            print 'Trying {} combination(s) of sites'.format(n_combos)
            combos = itertools.product(*current)

            useable = []
            time_init = time.time()
            for i, combo in enumerate(combos):
                if (time.time() - time_init) > 10:
                    print '{:.3f} percent complete'.format(float(i) / n_combos)
                    time_init = time.time()
                if check_useable(combo):
                    useable.append(combo)
            if useable:
                break

            position += 1

    if not exhaustive:
        current = [[x[0]] for x in sorted_walked]
        n_combos = 1
        for combo in current:
            n_combos *= len(combo)
        if not spannable(current):
            raise Exception('Couldn\'t do exhaustive search')

    sums = [sum([y[2] for y in x]) for x in useable]
    best = sums.index(max(sums))
    return useable[best]


def trim_directionally(pre_combo_list, max_len, direction):
    '''
    Traverse list of lists of (start, end, score) tuples of the type
    returned by StructureWindows. Removes any entries that cannot bridge the
    gap to the next set of overlaps (distance greater than max_distance).
    This greatly speeds up computation time by removing impossible combinations
    before they are evaluated computationally.

    :param pre_combo_list: list of lists of the output to StructureWindows,
                           representing a collection of each potential
                           overlap's evaluated windows.
    :type pre_combo_list: list
    :param max_len: Maximum length of each split-up output sequence.
    :type max_len: int

    '''

    pcl = [x for x in pre_combo_list]
    pcl_range = range(len(pcl) - 1)
    starts_list_list = [[y[0] for y in x] for x in pcl]
    for starts_list in starts_list_list:
        starts_list.sort()
    stops_list_list = [[y[1] for y in x] for x in pcl]
    for stops_list in stops_list_list:
        stops_list.sort()

    if direction == 'left':
        pcl_range.reverse()
    for i in pcl_range:
        starts_list = starts_list_list[i]
        stops_list = stops_list_list[i + 1]
        if direction == 'right':
            for j, start in enumerate(starts_list):
                spanned = [(x - start) <= max_len for x in stops_list]
                if not any(spanned):
                    pcl[i].pop(j)
                    return True
        elif direction == 'left':
            for j, stop in enumerate(stops_list):
                spanned = [(stop - x) <= max_len for x in starts_list]
                if not any(spanned):
                    pcl[i + 1].pop(j)
                    return True
    return False


def remove_nonspanning(pre_combo_list, max_len):
    '''
    Trims list of lists of (start, end, score) tuples to those
    that actually have a chance of spanning the sequence

    :param pre_combo_list: list of lists of the output to StructureWindows,
                           representing a collection of each potential
                           overlap's evaluated windows.
    :type pre_combo_list: list
    :param max_len: Maximum length of each split-up output sequence.
    :type max_len: int

    '''

    count = 1
    while True:
        left = trim_directionally(pre_combo_list, max_len=max_len,
                                  direction='right')
        right = trim_directionally(pre_combo_list, max_len=max_len,
                                   direction='left')
        if not left and not right:
            break

        count += 1
