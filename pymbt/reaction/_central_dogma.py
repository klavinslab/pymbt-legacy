'''
Reactions for the central dogma of biology - transcription and translation.

'''


from pymbt import sequence
from pymbt.reaction import utils


class Transcription(object):
    '''
    Transcription reaction - converts DNA to (m)RNA. Does not support
    post-transcriptional processing.

    '''

    def __init__(self, dna):
        '''
        :param dna: DNA input sequence.
        :type dna: DNA

        '''

        self.dna = dna
        self.rna = None  # yielded by .run()
        self.coding_rna = None  # yielded by .get_coding_rna()

    def run(self):
        '''
        Run the reaction.

        '''

        self.rna = utils.convert_sequence(self.dna, 'dna', 'rna')

        return self.rna

    def get_coding_rna(self):
        if not self.rna:
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
            if not valid[0]:
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

        if not valid[0]:
            raise Exception('Sequence has no start codon.')
        elif not stop:
            raise Exception('Sequence has no stop codon.')
        else:
            self.coding_rna = self.rna[start:stop]

        return self.coding_rna


class Translation(object):
    '''
    Translation reaction - converts (m)RNA to peptide. Does not support
    post-transcriptional/post-translational processing.

    '''

    def __init__(self, rna):
        '''
        :param dna: RNA input sequence.
        :type dna: RNA

        '''

        self.rna = rna
        self.coding_rna = None  # yielded by .get_coding_rna()
        self.peptide = None  # yielded by .run()
        self.coding_peptide = None  # yielded by .get_coding_peptide()

    def run(self):
        '''
        Run the reaction.

        '''

        self.peptide = utils.convert_sequence(self.rna,
                                              'rna',
                                              'peptide')

        return self.peptide

    # REFACTOR: this is totally redundant with the Transcription equivalent
    def get_coding_rna(self):
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
            if not valid[0]:
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

        if not valid[0]:
            raise Exception('Sequence has no start codon.')
        elif not stop:
            raise Exception('Sequence has no stop codon.')
        else:
            self.coding_rna = self.rna[start:stop]

        return self.coding_rna

    def get_coding_peptide(self):
        if not self.coding_rna:
            self.get_coding_rna()
        self.coding_peptide = utils.convert_sequence(self.coding_rna,
                                                     'rna',
                                                     'peptide')

        return self.coding_peptide


class ReverseTranscription(object):
    '''
    Reverse transcription reaction - converts (m)RNA to DNA.

    '''

    def __init__(self, rna):
        '''
        :param dna: RNA input sequence.
        :type dna: RNA

        '''

        self.rna = rna
        self.dna = None  # yielded by .run()

    def run(self):
        '''
        Run the reaction.

        '''

        self.dna = utils.convert_sequence(self.rna, 'rna', 'dna')

        return self.dna
