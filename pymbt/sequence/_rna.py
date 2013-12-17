'''RNA module.'''
import pymbt.reaction
from pymbt.sequence._sequence import BaseSequence
from pymbt.sequence import utils


class RNA(BaseSequence):
    '''RNA sequence.'''
    def __init__(self, rna, run_checks=True):
        '''
        :param rna: Input sequence (RNA).
        :type rna: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: pymbt.sequence.RNA instance.

        '''
        super(RNA, self).__init__(rna, 'rna', run_checks=run_checks)

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely-editable copy of the current sequence.
        :rtype: pymbt.sequence.RNA

        '''
        # Significant performance improvements by skipping alphabet check
        return type(self)(self._sequence, run_checks=False)

    def reverse_complement(self):
        '''Reverse complement sequence.

        :returns: The reverse-complement of the current sequence.
        :rtype: pymbt.sequence.RNA

        '''
        new_instance = self.copy()
        new_instance._sequence = utils.reverse_complement(self._sequence,
                                                          'rna')
        return new_instance

    def reverse_transcribe(self):
        '''Reverse transcribe to DNA.

        :returns: The reverse transcribed (DNA) version of the current RNA.
        :rtype: pymbt.sequence.DNA

        '''
        return pymbt.reaction.reverse_transcribe(self)

    def translate(self):
        '''Translate sequence into a peptide.

        :returns: A translated peptide from the current sequence.
        :rtype: pymbt.sequence.Peptide

        '''
        return pymbt.reaction.translate(self)

    def __contains__(self, query):
        """Defines `query in sequence` operator.

        :param query: query string or DNA sequence
        :type query: str or pymbt.sequence.DNA

        """
        return super(RNA, self).__contains__(query, "n")
