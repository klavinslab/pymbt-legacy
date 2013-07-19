'''The central dogma of biology - transcription and translation.'''
from pymbt import sequence
from pymbt.reaction import utils


def transcribe(dna):
    '''Transcribe DNA to RNA (no post-transcriptional processing).

    :param seq: Sequence to transcribe (DNA).
    :type seq: pymbt.sequence.DNA

    '''
    return utils.convert_sequence(dna, 'rna')


def translate(rna):
    '''Translate RNA to peptide.

    :param rna: Sequence to translate (RNA).
    :type rna: pymbt.sequence.RNA

    '''
    return utils.convert_sequence(rna, 'peptide')


def reverse_transcribe(rna):
    '''Reverse transcribe RNA to DNA.

    :param rna: Sequence to reverse transcribe (RNA).
    :type rna: pymbt.sequence.RNA

    '''
    return utils.convert_sequence(rna, 'dna')


def coding_sequence(rna):
    '''Extract coding sequence from an RNA template.

    :param seq: Sequence from which to extract a coding sequence.
    :type seq: pymbt.sequence.RNA
    :param material: Type of sequence ('dna' or 'rna')
    :type material: str

    '''
    if isinstance(rna, sequence.DNA):
        rna = transcribe(rna)
    codons_left = len(rna) // 3
    start_codon = sequence.RNA('aug')
    stop_codons = [sequence.RNA('uag'), sequence.RNA('uga'),
                   sequence.RNA('uaa')]
    start = None
    stop = None
    valid = [None, None]
    index = 0
    while codons_left:
        codon = rna[index:index+3]
        if valid[0] is None:
            if codon in start_codon:
                start = index
                valid[0] = True
        else:
            if codon in stop_codons:
                stop = index + 3
                valid[1] = True
                break
        index += 3
        codons_left -= 1

    if valid[0] is None:
        raise Exception('Sequence has no start codon.')
    elif stop is None:
        raise Exception('Sequence has no stop codon.')
    coding_rna = rna[start:stop]

    return coding_rna
