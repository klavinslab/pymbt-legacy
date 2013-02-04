import os # readres

from Bio import SeqIO # readref, readres

from pymbt.alignment.needle import needle

def readref(filepath):
    '''Helper function for reading reference genbank file'''
    sequence = SeqIO.read(path,'genbank')
    return(sequence)

#TODO: make this read ab1 files, and prefer them or something
def readres(dirpath):
    '''Helper function for reading .seq results files'''
    paths = [ x for x in os.listdir(dirpath) if x.endswith('.seq') ]
    sequences = [ SeqIO.read(x) for x in paths ]
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
        self.alignments = [ needle(ref,x) for x in res ]

        # TODO: score alignments. If score is too low, try rev-complement

        # Extract aligned sequences
        a_ref = self.alignments[0][0].seq.tostring()
        a_res = [ x[1].seq.tostring() for x in self.alignments ]
        self.aligned = {'reference':a_ref,'results':a_res}

        self.gaps = ''.join([ findgap([ y[i] for y in a_res ]) for i,x in enumerate(a_ref) ])

    def report():
        pass

    def plot():
        pass

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


