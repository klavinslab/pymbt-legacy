'''The central dogma of biology - transcription and translation.'''
from pymbt import sequence
from pymbt.reaction import utils


class Transcription(object):
    '''Transcribe DNA to RNA (no post-transcriptional processing).'''
    def __init__(self, dna):
        '''
        :param seq: Sequence to transcribe (DNA).
        :type seq: pymbt.sequence.DNA

        '''
        self.dna = dna
        self.rna = None

    def run(self):
        '''Transcribe the DNA sequence.'''
        self.rna = utils.convert_sequence(self.dna, 'rna')
        return self.rna

    def get_coding_rna(self):
        '''Extract coding RNA from sequence.'''
        if self.rna is None:
            self.run()
        codons_left = len(self.rna) // 3
        start_codon = sequence.RNA('aug')
        stop_codons = [sequence.RNA('uag'), sequence.RNA('uga'),
                       sequence.RNA('uaa')]
        start = None
        stop = None
        valid = [None, None]
        index = 0
        # REFACTOR: this code is super redundant but functions exactly right
        while codons_left:
            codon = self.rna[index:index+3]
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
        self.coding_rna = self.rna[start:stop]

        return self.coding_rna


class Translation(object):
    '''Translate RNA to peptide.'''
    def __init__(self, rna):
        '''
        :param rna: Sequence to translate (RNA).
        :type rna: pymbt.sequence.RNA

        '''
        self.rna = rna
        self.coding_rna = None
        self.peptide = None
        self.coding_peptide = None

    def run(self):
        '''Translate.'''
        self.peptide = utils.convert_sequence(self.rna, 'peptide')
        return self.peptide

    # REFACTOR: this is totally redundant with the Transcription equivalent
    def get_coding_rna(self):
        '''Extract coding RNA from sequence.'''
        codons_left = len(self.rna) // 3
        start_codon = sequence.RNA('aug')
        stop_codons = [sequence.RNA('uag'), sequence.RNA('uga'),
                       sequence.RNA('uaa')]
        start = None
        stop = None
        valid = [None, None]
        index = 0
        # REFACTOR: this code is super redundant but functions exactly right
        while codons_left:
            codon = self.rna[index:index+3]
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
        else:
            self.coding_rna = self.rna[start:stop]

        return self.coding_rna

    def get_coding_peptide(self):
        '''Extract coding peptide from sequence.'''
        if self.coding_rna is None:
            self.get_coding_rna()
        self.coding_peptide = utils.convert_sequence(self.coding_rna,
                                                     'peptide')

        return self.coding_peptide


class ReverseTranscription(object):
    '''Reverse transcribe RNA to DNA.'''
    def __init__(self, rna):
        '''
        :param rna: Sequence to reverse transcribe (RNA).
        :type rna: pymbt.sequence.RNA

        '''
        self.rna = rna
        self.dna = None  # yielded by .run()

    def run(self):
        '''Reverse transcribe.'''
        self.dna = utils.convert_sequence(self.rna, 'dna')
        return self.dna
