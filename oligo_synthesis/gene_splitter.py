import math
import itertools
import time

from pymbt.oligo_synthesis.structure_windows import context_walk

# BIG PROBLEM:
# context_walk results at some point are different than calling it directly
# TODO: circular (e.g. plasmid) version.
#   idea: find several good places to set as the origin around the plasmid,
#   then run split_gene using 'linearized' plasmid with that
#   60-base origin flanking either side. keep the best one.
#   this will probably take ~3-4 times as long as a linear fragment
# TODO: 200 minimum window of potential overlap was chosen arbitrarily
#   should get a sense of how good/bad the score is, *then* decide


def split_gene(seq,
               max_len=1100,
               min_context=200,
               core=60,
               context=90,
               step=10,
               force_exhaustive=False):

    n = 1
    while True:
        if len(seq) < n * max_len - (n - 1) * core:
            break
        else:
            n += 1
#    n = int(math.ceil(len(seq) / float(max_len)))

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
            stops = [x if x < len(seq) else len(seq) for x in stops]
            starts = [len(seq) - x for x in stops]
            starts = [x if x > 0 else 0 for x in starts]
            starts.reverse()
            olap_start_stop = [[starts[i], stops[i]] for i in range(len(stops))]
            olaps = [seq[x[0]:x[1]] for x in olap_start_stop]
            olap_w_context_start_stop = [[max(x[0]-context,0), min(x[1]+context,len(seq))] for x in olap_start_stop]
            olaps_w_context = [seq[x[0]:x[1]] for x in olap_w_context_start_stop]
            olap_w_context_lens = [x[1]-x[0] for x in olap_w_context_start_stop]
            if not any([x < min_context for x in olap_w_context_lens]):
                break 
            n += 1
    else:
        raise ValueError('Number of pieces is somehow zero.\
                         That shouldn\'t happen.')

    # TODO: see whether 'context' is screwing things up here - 
    # it doesn't need to be potential overlap
    walked_raw = []
    for i, x in enumerate(olaps_w_context):
        print('Analyzing %i of %i overlap(s).' % (i + 1, len(olaps)))
        walked_raw.append(context_walk(x, core, context, step, report=True))

    walked = []
    for i, x in enumerate(walked_raw):
        temp_list = []
        for y in x:
            y0 = y[0] + olap_w_context_start_stop[i][0]
            y1 = y[1] + olap_w_context_start_stop[i][0]
            y2 = y[2]
            temp_list.append((y0, y1, y2))
        walked.append(temp_list)

    score_start = [[starts[i] + y[0] for y in x] for i, x in enumerate(walked)]
    score_stop = [[stops[i] + y[1] for y in x] for i, x in enumerate(walked)]
    scores = [[y[2] for y in x] for x in walked]

    summary = {'starts': score_start,
               'stops': score_stop,
               'scores': scores}

    best_overlaps = find_best(walked, max_len, len(seq))
    olap_starts = [x[0] for x in best_overlaps]
    olap_stops = [x[1] for x in best_overlaps]
    scores = [x[2] for x in best_overlaps]
    olap_start_stop = [(x, olap_stops[i]) for i, x in enumerate(olap_starts)]

    best_starts = [0] + [x[0] for x in best_overlaps]
    best_stops = [x[1] for x in best_overlaps] + [len(seq)]
    b1 = best_starts
    b2 = best_stops
    final_seqs = [seq[b1[i]:b2[i]] for i in range(len(b1))]

    # TODO: return more than just sequences - at least start/stop of overlaps

    return {'sequences':final_seqs, 'overlaps':olap_start_stop, 'scores':scores}


def find_best(walked, max_distance, seq_len, force_exhaustive=False):
    # TODO BREAKTHROUGH:
    # new strategy:
    #   1. Similar strategy to before - increase possible sites number until
    #   sequence is spannable.
    #   2. Reduce the space:
    #       a. iterate through possibilities - exclude ones that will strictly
    #       never span. On edge this is simplest to understand - if a site in
    #       overlap 1 can't reach any of the sites in overlap 2, remove it.
    #       Still possible in the middle - can it span left-right?
    #       b. Use information acquired in last iteration?
    #       It wasn't possible to span before, but is by using
    #       at least one of the new sites. Is one more responsible for this than the others?
    #       c. Maybe go back a step, just add one at a time and see if it spans

    # TODO: when adding more positions to check for spanning,
    # add them intelligently instead of all across the board
    # alternatively, store all spanning options during a pass,
    # then pick the best one (probably a better idea)
    # For now, it goes with the first one it finds (not quite arbitrary)
    # For large fragments, this is especially important - 
    # e.g. Rob's 10kb fragment (20130225) would take ~3.5 hours
    # with 8 cores
    # to analyze using all combinations that have the potential
    # to minimally span the region

    # Use output of spannable() to add new entry at non-spannable site 
    # rather than increasing m incrementer
    # That way, a region without a connection gets more entries more quickly
    # and fewer combinations need to be calculated
    # Make it deterministic, though

    # Need to decide how to choose the one to which to add entries
    # Ideally, should perhaps be the one with the most room to work with
    # - the one with the least distance between it and its other, spannable
    # neighbor. If its neighbor isn't spannable, then increment it (bias to the left at first)

    # TODO:
    # should use exhaustive search when number of combinations is low,
    # above method when high.
    # 10 million cutoff?

    # TODO:
    # enable multiprocessing

    # Input is output of 'walked' adjusted for absolute sequence position
    # Output is indices of overlaps + score
    # (same format as entry of walked list)

    # Find the best combination of scores overlap positions given constraint:
    #   Synthesized pieces to be combined in Gibson must be < max_distance

    # strategy:
    # pick top score from each list, find all combinations (starts with 1)
    # see whether combination(s) have less than max_distance between positions
    # if not, add the next top score to the possibilites, and repeat last step

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
        for i in range(len(starts)-1):
            working_stops = [x for x in stops[i] if ((x - current_start) <= max_distance)]
            next_starts = [starts[i][j] for j, x in enumerate(stops[i]) if ((x - current_start) <= max_distance)]
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
    for i in range(3):
        pass
        #print([x[0:(i+1)] for x in sorted_walked])
    while exhaustive:
        current = [x[0:m] for x in sorted_walked]
        if not spannable(current):
            m += 1
        else:
            print([len(x) for x in current])
            remove_nonspanning(current, max_distance)
            print([len(x) for x in current])
            n_combos = 1
            for x in current:
                n_combos *= len(x)
            print(n_combos)
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
        if not spannable (current):
            # find first non-spanned position
            # try from the other side as well
            # if there is only one non-spanned position, add another entry for the one closest to its spanned neighbor
            # else, add another entry to the first non-spanned position
            # TODO: the above strategy is not guaranteed to work
            # alternative: alternate between above strategy and adding an entry to the right of the shortest-distance
            # spanner
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
        left = remove_from_left(pre_combo_list, max_len=max_len) # edits list in-place
        right = remove_from_right(pre_combo_list, max_len=max_len) # edits list in-place
        if not left and not right:
            both_false = True
