import csv
import itertools
import math
import time

from pymbt.oligo_synthesis.structure_windows import context_walk

# TODO: circular (e.g. plasmid) version.
#   1. Convert plasmid to 'linear' version - add 60bp from front to end
#   2. Evaluate as if linear
#   3. Rotate sequence by the window size and repeat until whole plasmid
#   has been analyzed
#   4. Pick the best one
#   Alternative:
#       Since we'd be analyzing the whole sequences anyways, just do that
#       right off the bat, store it, and try making spanning combinations


class GeneSplitter:
    '''A class that splits a large (~<10kb) sequence into smaller
    ones that are easier to clone, either via oligo assembly or PCR'''

    def __init__(self,
                 seq,
                 max_len=1100,
                 min_context=200,
                 core=60,
                 context=90,
                 step=10,
                 force_exhaustive=False):

        split = split_gene(seq,
                           max_len=max_len,
                           min_context=min_context,
                           core=core,
                           context=context,
                           step=step,
                           force_exhaustive=force_exhaustive)
        self.sequences = split['sequences']
        self.scores = split['scores']
        self.overlaps = split['overlaps']

    def write(self, path):
        f = open(path, 'w')
        cr = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
        cr.writerow(['sequence', 'olap_start', 'olap_stop', 'olap_score'])
        sta = [x[0] for x in self.overlaps] + ['NA']
        sto = [x[1] for x in self.overlaps] + ['NA']
        for i, x in enumerate(self.sequences):
                cr.writerow([x, sta[i], sto[i], self.scores[i]])


def split_gene(seq,
               max_len=1100,
               min_context=200,
               core=60,
               context=90,
               step=10,
               force_exhaustive=False):

    seq_len = len(seq)
    n = 1
    while True:
        if seq_len < n * max_len - (n - 1) * core:
            break
        else:
            n += 1

    # calculate potential overlap regions
    if n == 1:
        # doesn't need fancy calculations - no overlaps
        return seq
    elif n > 1:
        while True:
            # Generate overlap indices.
            # If overlaps smaller than min_context, increment n
            # Incrementing n in this way can lead to huge areas
            # for potential overlaps - may want to trim down for
            # redundant nupack usage
            stops = []
            for i in range(n - 1):
                stops.append((i + 1) * max_len - i * core)
            stops = [x if x < seq_len else seq_len for x in stops]
            starts = [seq_len - x for x in stops]
            starts = [x if x > 0 else 0 for x in starts]
            starts.reverse()
            olap_start_stop = [[starts[i], x] for i, x in enumerate(stops)]
            olaps = [seq[x[0]:x[1]] for x in olap_start_stop]
            with_context = []
            for x in olap_start_stop:
                start_i = max(x[0] - context, 0)
                stop_i = min(x[1] + context, seq_len)
                with_context.append([start_i, stop_i, seq_len])
            olaps_w_context = [seq[x[0]:x[1]] for x in with_context]
            olap_w_context_lens = [x[1] - x[0] for x in with_context]
            if not any([x < min_context for x in olap_w_context_lens]):
                break
            n += 1
    else:
        raise ValueError('Number of pieces is somehow zero.\
                         That shouldn\'t happen.')

    walked_raw = []
    for i, x in enumerate(olaps_w_context):
        print('Analyzing %i of %i overlap(s).' % (i + 1, len(olaps)))
        walked_raw.append(context_walk(x, core, context, step, report=True))

    walked = []
    for i, x in enumerate(walked_raw):
        temp_list = []
        for y in x:
            y0 = y[0] + with_context[i][0]
            y1 = y[1] + with_context[i][0]
            y2 = y[2]
            temp_list.append((y0, y1, y2))
        walked.append(temp_list)

    score_start = [[starts[i] + y[0] for y in x] for i, x in enumerate(walked)]
    score_stop = [[stops[i] + y[1] for y in x] for i, x in enumerate(walked)]
    scores = [[y[2] for y in x] for x in walked]

    summary = {'starts': score_start,
               'stops': score_stop,
               'scores': scores}

    best_overlaps = find_best(walked, max_len, seq_len)
    olap_starts = [x[0] for x in best_overlaps]
    olap_stops = [x[1] for x in best_overlaps]
    scores = [x[2] for x in best_overlaps]
    olap_start_stop = [(x, olap_stops[i]) for i, x in enumerate(olap_starts)]

    best_starts = [0] + [x[0] for x in best_overlaps]
    best_stops = [x[1] for x in best_overlaps] + [seq_len]
    b1 = best_starts
    b2 = best_stops
    final_seqs = [seq[b1[i]:b2[i]] for i in range(len(b1))]

    return {'sequences': final_seqs,
            'overlaps': olap_start_stop,
            'scores': scores}


def find_best(walked, max_distance, seq_len, force_exhaustive=False):
    '''Input to find_best is output of 'walked' adjusted for absolute sequence
    position.
    Output is indices of overlaps + score (the same format as entry of walked
    list)
    '''

    # sort results by score
    def sort_tuple(tuple_in):
        return sorted(tuple_in, key=lambda score: score[2], reverse=True)

    def spannable(start_stop_list):
        starts = [[y[0] for y in x] for x in start_stop_list]
        stops = [[y[1] for y in x] for x in start_stop_list]
        for x in stops:
            x.sort()
        for x in starts:
            x.sort()

        current_start = starts[0][-1]
        for i in range(len(starts) - 1):
            working_stops = []
            for x in stops[i]:
                if (x - current_start) <= max_distance:
                    working_stops.append(x)
            next_starts = []
            for j, x in enumerate(stops[i]):
                if (x - current_start) <= max_distance:
                    next_starts.append(starts[i][j])
            if not working_stops:
                return False
            else:
                current_start = next_starts[-1]
        return True

    def check_useable(tuple_in):
        for i in range(len(tuple_in) - 1):
            if tuple_in[i + 1][1] - tuple_in[i][0] > max_distance:
                return False
        return True

    sorted_walked = [sort_tuple(x) for i, x in enumerate(walked)]

    exhaustive = True
    m = 1
    while exhaustive:
        current = [x[0:m] for x in sorted_walked]
        if not spannable(current):
            m += 1
        else:
            remove_nonspanning(current, max_distance)
            n_combos = 1
            for x in current:
                n_combos *= len(x)
            if n_combos > 10000000 and not force_exhaustive:
                exhaustive = False
                break
            print('Trying %s combinations of top-scoring sites') % n_combos
            combos = itertools.product(*current)

            useable = []
            time_init = time.time()
            for i, combo in enumerate(combos):
                if (time.time() - time_init) > 10:
                    print('%.3f percent complete' % (float(i) / n_combos))
                    time_init = time.time()
                if check_useable(combo):
                    useable.append(combo)
            if useable:
                break

            m += 1

    if not exhaustive:
        current = [[x[0]] for x in sorted_walked]
        n_combos = 1
        for x in current:
            n_combos *= len(x)
        if not spannable(current):
            raise Exception('Couldn\'t do exhaustive search')

    sums = [sum([y[2] for y in x]) for x in useable]
    best = sums.index(max(sums))
    return useable[best]


def remove_from_left(pre_combo_list, max_len):
    '''Traverses list of lists of (start, end, score) tuples of the type
    returned by context_walk. Removes any entries that cannot bridge the
    gap to the next set of overlaps (distance greater than max_distance)'''

    pcl = pre_combo_list
    pcl_new = [x for x in pcl]
    for i in range(len(pcl) - 1):
        starts = [x[0] for x in pcl[i]]
        ends = [x[1] for x in pcl[i + 1]]

        for j, s in enumerate(starts):
            spanned = [(x - s) <= max_len for x in ends]
            if not any(spanned):
                pcl_new[i].pop(j)
                return True
    return False


def remove_from_right(pre_combo_list, max_len):
    '''Does the same thing as remove_from_left but in the reverse direction.
    These two functions have to be paired repeatedly in order to fully trim
    the sequence (remove_nonspanning will do this)'''

    pcl = pre_combo_list
    pcl_new = [x for x in pcl]
    # generate all indices but the first
    irange = range(len(pcl))
    irange.pop(0)
    irange.reverse()
    for i in irange:
        ends = [x[1] for x in pcl[i]]
        starts = [x[0] for x in pcl[i - 1]]

        for j, s in enumerate(ends):
            spanned = [(s - x) <= 1100 for x in starts]
            if not any(spanned):
                pcl_new[i].pop(j)
                return True
    return False


def remove_nonspanning(pre_combo_list, max_len):
    '''Trims list of lists of (start, end, score) tuples to those
    that actually have a chance of spanning the sequence'''

    both_false = False
    while not both_false:
        # these functions edit lists in place
        # lists are not copied on assignment in python, can be
        # painful to copy.
        # both functions edit the list in-place (no return values)
        left = remove_from_left(pre_combo_list, max_len=max_len)
        right = remove_from_right(pre_combo_list, max_len=max_len)
        if not left and not right:
            both_false = True
