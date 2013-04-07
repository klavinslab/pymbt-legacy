import os
import math

from matplotlib import pyplot
from matplotlib import cm
from Bio import SeqIO

from pymbt.alignment.needle import needle
from pymbt.sequence_manipulation import reverse_complement


class Sanger:
    '''Align and analyze Sanger sequencing results.'''
    def __init__(self, ref, res):
        '''
        :param ref: Reference sequence.
        :type ref: str.
        :param res: Sequencing result string. A list of strings is also valid.
        :type res: str.

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
        for i, x in enumerate(self.scores):
            if x < 1000:
                new_needle = needle(ref, reverse_complement(res[i]))
                self.alignments[i] = new_needle[0]
                self.scores[i] = new_needle[1]
        self.needle = []
        for i, x in enumerate(self.alignments):
            self.needle.append((self.alignments[i], self.scores[i]))

        # TODO:
        # 3. consensus / master sequence for plotting / report / analysis

        # Extract aligned sequences
        # Actually, need n ref sequences in case of insertions
        # - ref would have - inserted in that case.
        a_refs = [x[0].seq.tostring().upper() for x in self.alignments]
        a_res = [x[1].seq.tostring().upper() for x in self.alignments]

        for i, x in enumerate(a_refs):
            frontcount = 0
            backcount = 0
            while True:
                if x[-backcount - 1] == '-':
                    backcount += 1
                else:
                    break
            while True:
                if x[frontcount] == '-':
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
        for i, x in enumerate(a_res):
            self.gaps.append(''.join([_findgap([z]) for j, z in enumerate(x)]))
        self.ranges = []
        for x in self.gaps:
            first = x.find('X')
            last = (len(x) - x[::-1].find('X'))
            self.ranges.append((first, last))

        #mismatches
        self.mismatches = {}
        self.insertions = {}
        self.deletions = {}
        for i, x in enumerate(a_refs):
            curkey = str(i)
            #curkey = str(i) + ':' + self.resnames[i]
            for j, y in enumerate(x):
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

        for key, value in self.deletions.iteritems():
            last = (-2, -2)
            newlist = []
            for j, y in enumerate(value):
                if y[0] == last[1] + 1:
                    newlist[-1] = (last[0], y[1])
                else:
                    newlist.append(y)
                last = newlist[-1]
            self.deletions[key] = newlist

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
                for i, x in enumerate(value):
                    print ''
                    _sequences_display([self.ref[index], self.res[index]], x)
                    print ''
        if len(self.insertions) > 0:
            print '----------------------'
            print '----- INSERTIONS -----'
            print '----------------------'
            for key, value in self.insertions.iteritems():
                index = int(key)
                print '  ' + self.resnames[index]
                for i, x in enumerate(value):
                    print ''
                    _sequences_display([self.ref[index], self.res[index]], x)
                    print ''
        if len(self.deletions) > 0:
            print '---------------------'
            print '----- DELETIONS -----'
            print '---------------------'
            for key, value in self.deletions.iteritems():
                index = int(key)
                print '  ' + self.resnames[index]
                for i, x in enumerate(value):
                    print ''
                    _sequences_display([self.ref[index], self.res[index]], x)
                    print ''

    def plot(self):
        '''Plot visualization of results using matplotlib.'''
        # Step 1: turn results into ranges, bin those ranges for displaying
        ga = self.gaps
        bins = _disjoint_bins(self.ranges)

        # Step 2: plot 'reference' bar
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        gr0 = (0, len(self.aligned[0][0]))
        ax.broken_barh([gr0],
                       (13, 3),
                       facecolors='black',
                       edgecolors='none')

        # Plot and color features
        max_len = len(self.ref_raw.features)
        for i, x in enumerate(self.ref_raw.features):
            print i
            try:
                qual = x.qualifiers['label']
            except:
                qual = ''
            mid = int(math.ceil((x.location.start + x.location.end) / 2))
            locations = (x.location.start, x.location.end, mid, qual)
            ax.broken_barh([(locations[0], locations[1])],
                           (10, 9),
                           facecolors=cm.Set3(float(i) / max_len),
                           edgecolors='black')
            ax.text(locations[2] + 40,
                    15,
                    locations[3][0],
                    rotation=90)

        def add_discrepancies(index, bin):
            ax.plot(1000, 25)
            index = str(index)
            height = 22 + bin * 10

            for key, value in self.insertions.iteritems():
                if key == index:
                    for i, x in enumerate(value):
                        ax.plot(x[0], height, marker='o', color='k')
            for key, value in self.deletions.iteritems():
                if key == index:
                    for i, x in enumerate(value):
                        ax.plot(x[0], height, marker='^', color='k')
            for key, value in self.mismatches.iteritems():
                if key == index:
                    for i, x in enumerate(value):
                        ax.plot(x[0], height, marker='*', color='k')

        def wrap_name(str_in):
            n = 10
            out = [str_in[i:i + n] for i in range(0, len(str_in), n)]
            str_out = '\n'.join(out)
            return str_out

        # Step 3: plot results ranges
        for i, x in enumerate(bins):
            for y in x:
                index = y[2]
                gai = ga[index]
                gai_ind = gai.find('X')
                xy = (gai_ind, (len(gai) - gai[::-1].find('X')) - gai_ind)
                ax.broken_barh([xy],
                               (i * 10 + 20, 9),
                               facecolors='pink',
                               edgecolors='black')
                ax.text(xy[0] + 10,
                        i * 10 + 28,
                        wrap_name(self.resnames[y[2]]),
                        verticalalignment='top')
                add_discrepancies(index, i)

        #ax.set_ylim(5,35)
        ax.set_xlim(0, gr0[1])
        ax.set_xlabel('Base pairs from origin')
        ax.set_yticks([15, 25])
        ax.set_yticklabels(['Reference', 'Results'])
        ax.grid(True)
        ax.xaxis.grid(False)
        ax.yaxis.grid(False)

        pyplot.title('Alignment gap summary', fontstyle='normal')
        pyplot.show()

    def write_alignment():
        # TODO: custom format or standard (e.g. FASTA)? Implement both?
        pass

    def write_plot():
        pass

    def fixmismatch(result, position, newvalue):
        pass

    def fixindel():
        pass


def readref(filepath, ftype='genbank'):
    '''
    Read reference genbank file.

    :param filepath: Path to the file (typically genbank).
    :type filepath: str.
    :param ftype: Valid filetype readable by Bio.SeqIO.
    :type ftype: str.

    '''
    sequence = SeqIO.read(filepath, ftype)
    return sequence


def readres(dirpath):
    '''
    Read .seq results files from a dir.

    :param dirpath: Path to directory containing sequencing files.
    :type dirtpath: str.

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
    :type list_in: list.

    '''
    for x in list_in:
        if x != '-':
            return 'X'
    return '-'


def _sequences_display(seqs, start_stop, context=10, indent=4):
    '''
    Given two sequences to compare, display them and visualize non-matching
    regions.

    :param seqs: Sequences to display.
    :type seqs: list.
    :param start_stop: Indices to display.
    :type start_stop: tuple.
    :param context: Extra context to add on either side of the displayed
    sequences.
    :type context: int.
    :param indent: Indentation of the displayed text.
    :type indent: int.

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
    middle = ['|' if top[i] == bottom[i] else ' ' for i, x in enumerate(top)]
    middle = ''.join(middle)

    indent = ''.join([' ' for x in range(indent)])
    if start != stop:
        print indent + 'Positions %i to %i:' % (start, stop)
    else:
        print indent + 'Position %i:' % start
    print indent + top
    print indent + middle
    print indent + bottom


def _disjoint_bins(range_tuple_list):
    '''
    Construct disjoint bins given a list of range tuples.

    :param range_tuple_list: A list of tuples containing range values.
    :type range_tuple_list: list.

    '''
    rtl = range_tuple_list
    # number the ranges (third value in tuple)
    rtl = [(x[0], x[1], i) for i, x in enumerate(rtl)]

    # sort by range start
    rtl = sorted(rtl, key=lambda starts: starts[0])
    remaining = rtl

    done_binning = False
    n = 1
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
            n += 1

    return binned
