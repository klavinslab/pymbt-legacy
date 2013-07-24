'''Peptide module.'''
from pymbt.sequence._sequence import BaseSequence


class Peptide(BaseSequence):
    '''Peptide sequence.'''
    def __init__(self, peptide, run_checks=True):
        '''
        :param peptide: Input sequence (peptide).
        :type peptide: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool

        '''
        super(Peptide, self).__init__(peptide, 'peptide',
                                      run_checks=run_checks)

    def copy(self):
        '''Create a copy of the current instance.'''
        # Significant performance improvements by skipping alphabet check
        return type(self)(self._sequence, run_checks=False)
