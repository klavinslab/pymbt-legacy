# Functions to do basic DNA manipulations:
#   Reverse Complement
#   TODO: put more things here?
from string import maketrans, translate


# Simple reverse complement function
def reverse_complement(sequence):
    # TODO: verify that it's even DNA or RNA
    submat = maketrans('ATGCNatgcn', 'TACGNtacgn')
    sub_sequence = translate(sequence, submat)
    rev_sequence = sub_sequence[::-1]
    return(rev_sequence)
