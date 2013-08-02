'''Test sequencing module.'''
import os
from pymbt import analysis, seqio


class TestSanger(object):
    '''Test Sanger class.'''
    def __init__(self):
        ref_name = 'reference_sequence.gb'
        reference_path = os.path.join(os.path.dirname(__file__), ref_name)
        results_path = os.path.dirname(__file__)
        self.reference = self.read_reference(reference_path)
        self.results = self.read_results(results_path)

    def read_results(self, path):
        '''Read in sequencing results.'''
        seqs = seqio.read_sequencing(path)
        return seqs

    def read_reference(self, path):
        '''Read in sequencing reference sequence.'''
        seq = seqio.read_dna(path)
        return seq

    def test_sanger(self):
        '''Test that Sanger __init__ works without error.'''
        self.sanger = analysis.Sanger(self.reference, self.results)

    def test_alignment(self):
        '''Ensure that alignment is consistent.'''
        pass
