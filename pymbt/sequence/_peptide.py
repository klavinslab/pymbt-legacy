'''
Peptide object classes.

'''

import re
from pymbt.sequence import utils


class Peptide(object):
    '''
    Core Peptide sequence object.

    '''

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
        '''
        Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str

        '''

        pattern = pattern.lower()
        re_pattern = '(?=' + pattern + ')'
        indices_top = [index.start() for index in
                       re.finditer(re_pattern, self.peptide)]

        return indices_top

    def copy(self):
        '''
        Create a copy of the current instance.

        '''

        # Alphabet checking disabled on copy to improve performance
        new_instance = Peptide(self.peptide, run_checks=False)

        return new_instance

    def __getitem__(self, key):
        '''
        Indexing and slicing of sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object

        '''

        new_instance = self.copy()
        new_instance.peptide = self.peptide[key]

        return new_instance

    def __delitem__(self, index):
        '''
        Deletes sequence at index.

        param index: index to delete
        type index: int

        '''

        self.peptide = self.peptide[0:index] + self.peptide[index + 1:]

    def __setitem__(self, index, new_value):
        '''
        Sets index value to new value.

        '''

        insert = Peptide(new_value)
        new = self[0:index] + insert + self[index + 1:]
        self.peptide = new.peptide

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
        '''
        Coerce Peptide object to string. Only works for ungapped dsPeptide
        and top-strand ssPeptide.

        '''

        return self.peptide

    def __len__(self):
        '''
        Return length of all Peptide (including gaps) in object when built-in
        len function is used.

        '''

        return len(self.peptide)

    def __add__(self, other):
        '''
        Defines adding with + for Peptide objects.

        :param other: instance to be added to.
        :type other: Peptide

        '''

        tops = self.peptide + other.peptide
        new_instance = Peptide(tops, run_checks=False)

        return new_instance

    def __radd__(self, other):
        '''
        Add unlike types (enables sum function).

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
        else:
            # All is good - add the two DNA objects
            return self + other

    def __mul__(self, multiplier):
        '''
        Multiply Peptide by an integer to create concatenation.

        :param multiplier: Factor by which to multiply the sequence.
        :type multiplier: int

        '''

        # Input checking
        if multiplier != int(multiplier):
            msg = 'can\'t multiply sequence by non-integer.'
            raise TypeError(msg)

        # Test concatenation by adding once
        try:
            self + self
        except:
            raise Exception('Failed to add, so cannot multiply.')

        # Isolate sequence strands, multiply strings, recreate Peptide
        tops = self.peptide * multiplier
        new_instance = Peptide(tops, run_checks=False)

        return new_instance

    def __eq__(self, other):
        '''
        Test Peptide object equality.

        '''

        if vars(self) == vars(other):
            return True
        else:
            return False

    def __ne__(self, other):
        '''
        Test Peptide object equality.

        '''

        return not (self == other)
