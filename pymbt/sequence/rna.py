'''
RNA object classes.

'''

import math
import re
from pymbt.sequence import utils


class RNA(object):
    '''
    Core RNA sequence object.

    '''
    def __init__(self, sequence, run_checks=True):
        '''
        :param sequence: Input sequence (RNA).
        :type sequence: str
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool

        '''

        self.top = sequence
        if run_checks:
            self.top = utils.check_seq(sequence, 'rna')

        self.bottom = ''.join('-' for x in self.top)

    def reverse_complement(self):
        '''
        Reverse complement top and bottom strands.

        '''
        new_instance = self.copy()
        new_instance.top = utils.reverse_complement(self.top, 'rna')

        return new_instance

    def five_resect(self, n_bases=None):
        '''
        Remove bases from 5' end of top strand.


        :param n_bases: Number of bases cut back. Defaults to removing entire
                        top strand
        :type n_bases: int

        '''
        n_blanks = len(self.top) - len(self.top.lstrip('-'))
        if n_bases and n_bases < len(self.top):
            gap_length = n_blanks + n_bases
            new_gap = '-' * gap_length
            new_top = new_gap + self.top[gap_length:]
        else:
            new_top = '-' * len(self.top)

        new_instance = self.copy()
        new_instance.top = new_top
        new_instance.remove_end_gaps()

        return new_instance

    def three_resect(self, n_bases=None):
        '''
        Remove bases from 3' end of top strand.

        :param n_bases: Number of bases cut back. Defaults to removing entire
                        top strand
        :type n_bases: int

        '''
        n_blanks = len(self.top) - len(self.top.rstrip('-'))
        if n_bases and n_bases < len(self.top):
            gap_length = n_blanks + n_bases
            new_gap = '-' * gap_length
            new_top = self.top[:-gap_length] + new_gap
        else:
            new_top = '-' * len(self.top)

        new_instance = self.copy()
        new_instance.top = new_top
        new_instance.remove_end_gaps()

        return new_instance

    def locate(self, pattern):
        '''
        Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str

        '''

        pattern = pattern.lower()
        template_top = self.top
        template_bottom = self.bottom
        re_pattern = '(?=' + pattern + ')'
        indices_top = [index.start() for index in
                       re.finditer(re_pattern, template_top)]
        indices_bottom = [index.start() for index in
                          re.finditer(re_pattern, template_bottom)]
        # TODO: if pattern is inverted repeat, throw out top/bottom redundant
        # matches. For now will just check top strand only - but this will
        # fail if there's gaps

        def check_inv(pattern):
            '''
            Check whether pattern is palindrome.
            :param pattern: pattern to test.
            :type pattern: str

            '''
            p_len = len(pattern)
            wing = int(math.floor(p_len / 2))
            if p_len % 2 != 0:
                l_wing = pattern[0:wing + 1]
                r_wing = pattern[wing:]
            else:
                l_wing = pattern[0: wing]
                r_wing = pattern[wing:]
            if l_wing == utils.reverse_complement(r_wing, 'rna'):
                return True
            else:
                return False

        inverted_repeat = check_inv(pattern)
        if inverted_repeat:
            # subtract all occurrences in top from bottom
            subtract = [len(self.top) - index - len(pattern) for index in
                        indices_top]
            indices_bottom = [x for x in subtract if x not in indices_bottom]

        return (indices_top, indices_bottom)

    def copy(self):
        '''
        Create a copy of the current instance.

        '''

        # Note: alphabet checking disabled on copy to improve preformance -
        # will bottleneck many workflows that rely on slicing otherwise
        new_instance = RNA(self.top, run_checks=False)
        return new_instance

    def remove_end_gaps(self):
        '''
        Removes gaps from ends of the sequence.

        '''
        top = self.top
        top_nogaps = len(top.lstrip('-')), len(top.rstrip('-'))
        self.top = top_nogaps
        self.bottom = ''.join(['-' for x in self.top])

    def __getitem__(self, key):
        '''
        Indexing and slicing of sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object

        '''
        # TODO: throw proper error when index is out of range
        new_instance = self.copy()
        if isinstance(key, int):
            new_instance.top = new_instance.top[key]
            new_instance.bottom = new_instance.bottom[::-1][key]

            return new_instance
        elif isinstance(key, slice):
            new_instance.top = new_instance.top.__getitem__(key)
            bottom_rev = new_instance.bottom[::-1]
            new_instance.bottom = bottom_rev.__getitem__(key)[::-1]
            return new_instance

    def __delitem__(self, index):
        '''
        Deletes sequence at index.

        param index: index to delete
        type index: int

        '''

        new = self[0:index] + self[index + 1:]
        self.top = new.top
        self.bottom = new.bottom

    def __setitem__(self, index, new_value):
        '''
        Sets index value to new value.

        '''

        insert = RNA(new_value)
        new = self[0:index] + insert + self[index + 1:]
        self.top = new.top
        self.bottom = new.bottom

    def __repr__(self):
        '''
        String to print when object is called directly.

        '''
        show_bases = 40
        bottom = self.bottom[::-1]
        if len(self.top) < 90:
            to_print = [self.top, bottom]
        else:
            top = ''.join([self.top[0:show_bases], ' ... ',
                           self.top[-show_bases:]])
            bottom = ''.join([bottom[0:show_bases], ' ... ',
                              bottom[-show_bases:]])
            to_print = [top, bottom]
        first_line = 'RNA:'
        to_print = [first_line] + to_print

        return '\n'.join(to_print)

    def __str__(self):
        '''
        Coerce RNA object to string. Only works for ungapped dsRNA
        and top-strand ssRNA.

        '''
        # TODO: implement bottom-only strand ssRNA as well?
        gaps = sum([1 for x in self.top if x == '-'])
        #gaps += sum([1 for x in self.bottom if x == '-'])
        if gaps:
            # note: should really implement this - getting an error from a
            # print line is stupid
            msg1 = 'No string coercion method for gapped '
            msg2 = 'sequences.'
            print msg1 + msg2
            return NotImplemented
        else:
            return self.top

    def __len__(self):
        '''
        Return length of all RNA (including gaps) in object when built-in
        len function is used.

        '''
        return len(self.top)

    def __add__(self, other):
        '''
        Defines adding with + for RNA objects. Only works for ungapped dsRNA
        and top-only ssRNA.

        :param other: instance to be added to.
        :type other: compatible sequence object (currently only RNA).

        '''

        # Check self for terminal gaps:
        if len(self):
            self_gaps = (self.top[-1] == '-', self.bottom[0] == '-')
        else:
            self_gaps = (False, False)
        if len(other):
            other_gaps = (other.top[0] == '-', other.bottom[-1] == '-')
        else:
            other_gaps = (False, False)

        # Only time this should fail is when there's a discontinuity. A
        # discontinuity appears when the first and second entries for the gaps
        # tuples are exactly the opposite of each other
        if self_gaps[0] != other_gaps[0] and self_gaps[1] != other_gaps[1]:
            # TODO: pretty sure this is the wrong way to throw an exception
            # here
            print 'Discontinuity at ends of RNA objects - can\'t add.'
            return NotImplemented

        tops = self.top + other.top
        bottoms = other.bottom + self.bottom

        new_instance = self.copy()
        new_instance.top = tops
        new_instance.bottom = bottoms

        return new_instance

    def __mul__(self, multiplier):
        '''
        Multiply RNA by an integer to create concatenation.

        :param multiply: Factor by which to multiply the sequence.
        :type multiply: int

        '''

        # Input checking
        if multiplier != int(multiplier):
            msg = 'can\'t multiply sequence by non-int.'
            raise TypeError(msg)

        # Try adding once as a test. This is slow for low values of b
        # but very fast for higher values
        try:
            self + self
        except:
            raise Exception('Failed to add, so cannot multiply.')

        # If addition check passes, just isolate top and bottom strands, do
        # multiplication to the strings (fast), and recreate RNA
        tops = self.top * multiplier

        return RNA(tops, run_checks=False)
