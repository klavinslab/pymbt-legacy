# Functions to do basic DNA manipulations

from string import maketrans,translate

# Simple reverse complement function
def reverse_complement(sequence):
    submat = maketrans('ATGCatgc','TACGtacg')
    sub_sequence = translate(sequence,submat)
    rev_sequence = sub_sequence[::-1]
    return(rev_sequence)
