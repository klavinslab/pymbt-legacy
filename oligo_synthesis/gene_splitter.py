import math
import itertools

from pymbt.oligo_synthesis.structure_windows import context_walk

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
               step=10):

    eval_len = core + context
    # Trim down sequence to the necessary overlaps
    n = int(math.ceil(len(seq) / float(max_len)))

    # calculate potential overlap regions
    if n == 1:
        # doesn't need fancy calculations - no overlaps
        return seq
    elif n > 1:
        # ensure at least min_context bp in each potential overlap region
        while max_len - (len(seq) - (n - 1) * max_len) < min_context:
            n += 1
        starts = [len(seq) - (x + 1) * max_len for x in range(n - 1)][::-1]
        stops = [(x + 1) * max_len for x in range(n - 1)]
        olaps = [seq[starts[i]:x] for i, x in enumerate(stops)]
    else:
        raise ValueError('Number of pieces is somehow zero.\
                         That shouldn\'t happen.')

    walked_raw = []
    for i, x in enumerate(olaps):
        print('Analyzing %i of %i overlap(s).' % (i + 1, len(olaps)))
        walked_raw.append(context_walk(x, core, context, step, report=True))

    walked = []
    for i, x in enumerate(walked_raw):
        temp_list = []
        for y in x:
            y0 = y[0] + starts[i]
            y1 = y[1] + starts[i]
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

    best_starts = [0] + [x[0] for x in best_overlaps]
    best_stops = [x[1] for x in best_overlaps] + [len(seq)]
    b1 = best_starts
    b2 = best_stops
    final_seqs = [seq[b1[i]:b2[i]] for i in range(len(b1))]

    return final_seqs


def find_best(walked, max_distance, seq_len):
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

    sorted_walked = [sort_tuple(x) for i, x in enumerate(walked)]

    def check_useable(tuple_in):
        starts_i = [0] + [y[0] for y in tuple_in]
        stops_i = [y[1] for y in tuple_in] + [seq_len]
        for j, y in enumerate(starts_i):
            if stops_i[j] - y > max_distance:
                return False
        return True

    m = 1
    combo_found = False
    while not combo_found:
        current = [x[0:m] for x in sorted_walked]
        combos = list(itertools.product(*current))

        useable = [check_useable(x) for x in combos]

        combo_found = any(useable)
        m += 1
    # Trim to those that span the gene with max_distance condition
    combos = [x for i, x in enumerate(combos) if useable[i]]
    sums = [sum([y[2] for y in x]) for x in combos]
    best = sums.index(max(sums))
    return combos[best]
