'''
RNA object classes.

'''

import re
from pymbt.sequence import utils


class RNA(object):
    '''
    Core RNA sequence object.

    '''

    def __init__(self, seq, run_checks=True):
        '''
        :param seq: Input sequence (RNA).
        :type seq: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool

        '''

        if run_checks:
            self.top = utils.process_seq(seq, 'rna')
        else:
            self.top = seq

        self.bottom = ''.join('-' for x in self.top)

    def reverse_complement(self):
        '''
        Reverse complement sequence.

        '''

        new_instance = self.copy()
        new_instance.top = utils.reverse_complement(self.top, 'rna')

        return new_instance

    def five_resect(self, n_bases):
        '''
        Remove bases from 5' end of top strand.


        :param n_bases: Number of bases cut back.
        :type n_bases: int

        '''

        new_instance = self[n_bases::]

        return new_instance

    def three_resect(self, n_bases):
        '''
        Remove bases from 3' end of top strand.

        :param n_bases: Number of bases cut back.
        :type n_bases: int

        '''

        new_instance = self[:-n_bases]

        return new_instance

    def locate(self, pattern):
        '''
        Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str

        '''

        pattern = pattern.lower()
        re_pattern = '(?=' + pattern + ')'
        indices_top = [index.start() for index in
                       re.finditer(re_pattern, self.top)]

        return indices_top

    def copy(self):
        '''
        Create a copy of the current instance.

        '''

        # Alphabet checking disabled on copy to improve performance
        new_instance = RNA(self.top, run_checks=False)

        return new_instance

    def __getitem__(self, key):
        '''
        Indexing and slicing of sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object

        '''

        new_instance = self.copy()
        new_instance.top = self.top[key]
        new_instance.bottom = self.bottom[::-1][key][::-1]

        return new_instance

    def __delitem__(self, index):
        '''
        Deletes sequence at index.

        param index: index to delete
        type index: int

        '''

        top_list = list(self.top)
        bottom_list = list(self.top)
        top_list.pop(index)
        bottom_list.pop(index)
        self.top = ''.join(top_list)
        self.bottom = ''.join(top_list)

    def __setitem__(self, index, new_value):
        '''
        Sets index value to new value.

        '''

        top_list = list(self.top)
        insert = RNA(new_value)
        top_list[index] = str(insert)
        self.top = ''.join(top_list)

    def __repr__(self):
        '''
        String to print when object is called directly.

        '''

        show = 40
        bottom = self.bottom[::-1]

        if len(self.top) < 90:
            top = self.top
            bottom = '-' * len(top)
        else:
            top = ''.join([self.top[0:show], ' ... ', self.top[-show:]])
            bottom = ''.join(['-' * show, ' ... ', '-' * show])

        first_line = 'RNA:'
        to_print = '\n'.join([first_line, top, bottom])

        return to_print

    def __str__(self):
        '''
        Coerce RNA object to string. Only works for ungapped dsRNA
        and top-strand ssRNA.

        '''

        return self.top

    def __len__(self):
        '''
        Return length of all RNA (including gaps) in object when built-in
        len function is used.

        '''

        return len(self.top)

    def __add__(self, other):
        '''
        Defines adding with + for RNA objects.

        :param other: instance to be added to.
        :type other: RNA

        '''

        tops = self.top + other.top
        new_instance = RNA(tops, run_checks=False)

        return new_instance

    def __radd__(self, other):
        '''
        Add unlike types (enables sum function).

        :param self: object of self type.
        :type self: RNA
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
        '''
        Multiply RNA by an integer to create concatenation.

        :param multiplier: Factor by which to multiply the sequence.
        :type multiplier: int

        '''

        # Input checking
        if multiplier != int(multiplier):
            msg = 'can\'t multiply sequence by non-integer.'
            raise TypeError(msg)

        # Isolate top and bottom strands, multiply strings, recreate RNA
        tops = self.top * multiplier
        new_instance = RNA(tops, run_checks=False)

        return new_instance

    def __eq__(self, other):
        '''
        Test RNA object equality.

        '''

        if vars(self) == vars(other):
            return True
        else:
            return False

    def __ne__(self, other):
        '''
        Test RNA object equality.

        '''

        return not (self == other)

    def __contains__(self, pattern):
        '''
        Specially-defined `x in y` behavior

        '''

        if str(pattern) in str(self):
            return True
        else:
            return False
