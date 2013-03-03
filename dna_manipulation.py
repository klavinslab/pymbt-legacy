from string import maketrans, translate


# Simple reverse complement function
def reverse_complement(sequence):
    # TODO: verify that it's even DNA or RNA
    submat = maketrans('ATGCNatgcn', 'TACGNtacgn')
    sub_sequence = translate(sequence, submat)
    rev_sequence = sub_sequence[::-1]
    return rev_sequence


def check_alphabet(sequence, material='dna'):
    '''Verifies that a given sequence is made only of DNA or RNA'''
    if material is 'dna':
        err_msg = 'DNA'
        alphabet = 'ATGCatgc'
    elif material is 'rna':
        err_msg = 'RNA'
        alphabet = 'AUGCaugc'
    for char in sequence:
        if char not in alphabet:
            raise ValueError('Sequence has a non-%s character') % err_msg
    return sequence
