'''Peptide object classes.'''

import re
from pymbt.sequence import utils


class Peptide(object):
    '''Peptide sequence.'''

    def __init__(self, seq, run_checks=True):
        '''
        :param seq: Input sequence (peptide).
        :type seq: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool

        '''
        if run_checks:
            self.peptide = utils.process_seq(seq, 'peptide')
        else:
            self.peptide = seq

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str

        '''
        pattern = pattern.lower()
        re_pattern = '(?=' + pattern + ')'
        indices_top = [index.start() for index in
                       re.finditer(re_pattern, self.peptide)]

        return indices_top

    def copy(self):
        '''Create a copy of the current instance.'''
        # Alphabet checking disabled on copy to improve performance
        new_instance = Peptide(self.peptide, run_checks=False)

        return new_instance

    def __getitem__(self, key):
        '''Index and slice sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object

        '''
        new_instance = self.copy()
        new_instance.peptide = self.peptide[key]

        return new_instance

    def __delitem__(self, index):
        '''Delete sequence at index.

        param index: index to delete
        type index: int

        '''
        peptide_list = list(self.peptide)
        peptide_list.pop(index)
        self.peptide = ''.join(peptide_list)

    def __setitem__(self, index, new_value):
        '''Set value at index to new value.

        '''
        insert = Peptide(new_value)
        peptide_list = list(self.peptide)
        peptide_list[index] = str(insert)
        self.peptide = ''.join(peptide_list)

    def __repr__(self):
        '''
        String to print when object is called directly.

        '''
        show = 40

        if len(self.peptide) < 90:
            peptide = self.peptide
        else:
            peptide = ''.join([self.peptide[0:show], ' ... ',
                              self.peptide[-show:]])

        first_line = 'Peptide:'
        to_print = '\n'.join([first_line, peptide])

        return to_print

    def __str__(self):
        '''Coerce Peptide object to string.'''
        return self.peptide

    def __len__(self):
        '''Return length of the peptide'''
        return len(self.peptide)

    def __add__(self, other):
        '''Defines addition for Peptide objects.

        :param other: instance to be added to.
        :type other: Peptide

        '''
        tops = self.peptide + other.peptide
        new_instance = Peptide(tops, run_checks=False)
        return new_instance

    def __radd__(self, other):
        '''Add unlike types (enables sum function).

        :param self: object of self type.
        :type self: Peptide
        :param other: object of any other type.
        :param other: anything

        '''
        if other == 0:
            # sum(list) adds to zero first, so ignore it
            return self
        elif type(self) != type(other):
            # Ensure types are the same
            msg = 'unsupported operand type(s) for +: {} and {}'.format(self,
                                                                        other)
            raise TypeError(msg)

    def __mul__(self, multiplier):
        '''Multiply Peptide by an integer to create concatenation.

        :param multiplier: Factor by which to multiply the sequence.
        :type multiplier: int

        '''
        # Input checking
        if multiplier != int(multiplier):
            msg = 'can\'t multiply sequence by non-integer.'
            raise TypeError(msg)

        # Isolate sequence strands, multiply strings, recreate Peptide
        tops = self.peptide * multiplier
        new_instance = Peptide(tops, run_checks=False)

        return new_instance

    def __eq__(self, other):
        '''Test Peptide object equality.

        :param other: other peptide
        :type other: pymbt.sequence.Peptide

        '''
        if vars(self) == vars(other):
            return True
        else:
            return False

    def __ne__(self, other):
        '''Test Peptide object inequality.

        :param other: other peptide
        :type other: pymbt.sequence.Peptide

        '''
        return not (self == other)

    def __contains__(self, pattern):
        '''Specially-defined `x in y` behavior

        :param pattern: pattern to find.
        :type pattern: str

        '''
        if str(pattern) in str(self):
            return True
        else:
            return False
