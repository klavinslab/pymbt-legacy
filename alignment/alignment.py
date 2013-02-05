import os # readres

from matplotlib import pyplot # plotting
from Bio import SeqIO # readref, readres

from pymbt.alignment.needle import needle
from pymbt.dna_manipulation import reverse_complement

def readref(filepath,ftype='genbank'):
    '''Helper function for reading reference genbank file'''
    sequence = SeqIO.read(filepath,ftype).seq.tostring()
    return(sequence)

#TODO: make this read ab1 files, and prefer them or something
def readres(dirpath):
    '''Helper function for reading .seq results files'''
    paths = [ x for x in os.listdir(dirpath) if x.endswith('.seq') ]
    sequences = [ SeqIO.read(dirpath+x,'fasta').seq.tostring() for x in paths ]
    return(sequences)

def findgap(list_in):
    for x in list_in:
        if x != '-':
            return('X')
    return('-')

class Sanger:
    '''Sanger sequencing alignment generation and analysis class'''
    def __init__(self,ref,res):
        # make results a list if there's just one
        if type(res) is str:
            res = [res]
        # Reduce results to largest unambiguous segment and align
        split = [ x.split('N') for x in res ]
        lengths = [ [ len(y) for y in x ] for x in split ]
        max_positions = [ x.index(max(x)) for x in lengths ]
        res = [ x[max_positions[i]] for i,x in enumerate(split) ]
        self.needle = [ needle(ref,x) for x in res ]
        self.alignments = [ x[0] for x in self.needle ]
        self.scores = [ x[1] for x in self.needle ]
        for i,x in enumerate(self.scores):
            if x < 1000:
                new_needle = needle(ref,reverse_complement(res[i]))
                self.alignments[i] = new_needle[0]
                self.scores[i] = new_needle[1]

        # TODO: score alignments. If score is too low, try rev-complement

        # Extract aligned sequences
        # Actually, need n ref sequences in case of insertions
        # - ref would have - inserted in that case.
        a_refs = [ x[0].seq.tostring() for x in self.alignments ]
        a_res = [ x[1].seq.tostring() for x in self.alignments ]
        #self.aligned = {'reference':a_refs,'results':a_res}

        self.aligned = [ (a_refs[i],a_res[i]) for i,x in enumerate(a_refs) ]
        for i,x in enumerate(a_refs):
            for j,y in enumerate(x):
                if y == '-':
        #            print('gap found in' + str(i) + ':' + str(y))
                    pass

        # Need to find a *consensus* before displaying anything
        # Need to deal with predicted insertions, then: coordinate change


        self.gaps = [ ''.join([ findgap([z]) for j,z in enumerate(x) ]) for i,x in enumerate(a_res) ]

    def report():
        pass

    def plot(self):
        # Step 1: turn results into ranges, bin those ranges for displaying
        ga = self.gaps
        res_ranges = [ (x.find('X') , (len(x)-x[::-1].find('X'))) for x in ga ]
        bins = disjoint_bins(res_ranges)

        # Step 2: plot 'reference' bar
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        gr0 = (0,len(self.aligned[0][0]))

        ax.broken_barh([ gr0 ] , (10, 9), facecolors='blue',edgecolors='none')

        # Step 3: plot results ranges
        for i,x in enumerate(bins):
            for y in x:
                index = y[2]
                gai = ga[index]
                xy = (gai.find('X') , (len(gai)-gai[::-1].find('X'))-gai.find('X'))
                ax.broken_barh([ xy ] , (i*10+20, 9), facecolors='red',edgecolors='none')

        #ax.set_ylim(5,35)
        ax.set_xlim(0,gr0[1])
        ax.set_xlabel('Base pairs from origin')
        ax.set_yticks([15,25])
        ax.set_yticklabels(['Reference', 'Results'])
        ax.grid(True)
        
        pyplot.title('Alignment gap summary',fontstyle='normal')
        pyplot.show()

    def write_alignment():
        pass

    def write_plot():
        pass

    def mismatches():
        pass

    def indels():
        pass

    def fixmismatch():
        pass

    def fixindel():
        pass


def disjoint_bins(range_tuple_list):
    '''Construct disjoint bins given a list of range tuples'''
    rtl = range_tuple_list
    # number the ranges (third value in tuple)
    rtl = [ (x[0],x[1],i) for i,x in enumerate(rtl) ]

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
            nextbin += [ x for x in remaining if x[0] <= last[1] ]
            remaining = [ x for x in remaining if x[0] > last[1] ]
        binned.append(current_bin)
        if len(remaining) == 0 and len(nextbin) == 0:
            done_binning = True 
        else:
            remaining = nextbin
            n += 1

    return(binned)
    #take first range, then next closest that doesn't overlap, etc until done
    #put remainder in next bin
    #repeat until done
