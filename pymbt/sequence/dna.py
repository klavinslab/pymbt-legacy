import math
import re

from pymbt import sequence_utils


class DNA(object):
    def __init__(self, sequence, bottom=None, type='linear', stranded='ds',
                 features=[]):
        sequence_utils.check_alphabet(sequence)
        self.type = type
        self.top = sequence.lower()
        self.features = features
        # TODO: eliminate this attribute - just let it be implicit / checkable
        # with a method e.g. DNA.stranded()
        self.stranded = stranded
        self.name = ''

        if bottom:
            sequence_utils.check_alphabet(bottom)
            self.bottom = bottom
            # TODO: check for complementation between top/bottom if allowing
            # manual bottom strand input. Mismatched complexes should
            # not be implemented yet
        elif stranded == 'ss':
            self.bottom = ''.join('-' for x in self.top)
        elif stranded == 'ds':
            self.bottom = sequence_utils.reverse_complement(self.top)

    def reverse_complement(self):
        '''
        Reverse complement top and bottom strands.
        '''
        new_instance = self.copy()
        if self.stranded == 'ds':
            new_instance.top = self.bottom
            new_instance.bottom = self.top
        elif self.stranded == 'ss':
            new_instance.top = sequence_utils.reverse_complement(self.top)

        return new_instance

    def circularize(self):
        '''
        Circularize linear DNA.
        '''
        # Currently does nothing other than changing an attribute
        new_instance = self.copy()
        new_instance.type = 'circular'

        return new_instance

    def linearize(self, index=0):
        '''
        Linearize circular DNA at a specific index.
        '''
        # TODO: collections.deque makes this easier.
        if self.type == 'linear':
            raise Exception('Cannot relinearize linear DNA.')
        top_list = [base for base in self.top]
        bottom_list = [base for base in self.bottom]
        if index > 0:
            for x in range(index):
                top_pop = top_list.pop(0)
                top_list.append(top_pop)
                bottom_list = [bottom_list.pop()] + bottom_list
        elif index < 0:
            for x in range(abs(index)):
                top_list = [top_list.pop()] + top_list
                bottom_pop = bottom_list.pop(0)
                bottom_list.append(bottom_pop)

        new_instance = self.copy()
        new_instance.top = ''.join(top_list)
        new_instance.bottom = ''.join(bottom_list)
        new_instance.type = 'linear'

        return new_instance

    def five_resect(self, n_bases=None):
        '''
        Remove bases from 5' end of top strand.
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
        new_instance._remove_end_gaps()

        return new_instance

    def three_resect(self, n_bases=None):
        '''
        Remove bases from 3' end of top strand.
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
        new_instance._remove_end_gaps()

        return new_instance

    def locate(self, pattern):
        pattern = pattern.lower()
        if self.type == 'circular':
            # TODO: this probably doesn't account for gaps
            max_len = min(len(self.top), len(pattern)) + len(self.top) - 1
            template_top = self.top + self.top[0:max_len]
            template_bottom = self.bottom + self.bottom[0:max_len]
        else:
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
            p_len = len(pattern)
            wing = int(math.floor(p_len / 2))
            if p_len % 2 != 0:
                l_wing = pattern[0:wing + 1]
                r_wing = pattern[wing:]
            else:
                l_wing = pattern[0: wing]
                r_wing = pattern[wing:]
            if l_wing == sequence_utils.reverse_complement(r_wing):
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
        # Making a new instance this way is slow. It's just as slow if
        # initialized with an empty string '' due to alphabet checking
        new_instance = DNA(self.top, bottom=self.bottom, type=self.type,
                           stranded=self.stranded, features=self.features)
#        new_instance = copy.deepcopy(self)
        return new_instance

    def _remove_end_gaps(self):
        '''
        Removes double-stranded gaps from ends of the sequence.
        '''
        top = self.top
        bottom = self.bottom
        top_nogaps = len(top.lstrip('-')), len(top.rstrip('-'))
        bottom_nogaps = len(bottom.lstrip('-')), len(bottom.rstrip('-'))
        top_gaps = [len(top) - x for x in top_nogaps]
        bottom_gaps = [len(bottom) - x for x in bottom_nogaps]

        left_trim = min(top_gaps[0], bottom_gaps[1])
        right_trim = min(top_gaps[1], bottom_gaps[0])

        # TODO: should define slicing for DNA objects and use that here instead
        if left_trim:
            top = top[left_trim:]
            bottom = bottom[:-left_trim]
        if right_trim:
            top = top[:-right_trim]
            bottom = bottom[right_trim:]

        self.top = top
        self.bottom = bottom

    def __getitem__(self, key):
        '''
        Indexing and slicing of sequences.
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

        if new_instance.type == 'circular' and len(new_instance) != len(self):
            new_instance.type = 'linear'
            return new_instance

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
        first_line = '{} {}DNA:'.format(self.type, self.stranded)
        to_print = [first_line] + to_print

        return '\n'.join(to_print)

    def __str__(self):
        '''
        Coerce DNA object to string. Only works for ungapped dsDNA
        and top-strand ssDNA.
        '''
        # TODO: implement bottom-only strand ssDNA as well?
        gaps = sum([1 for x in self.top if x == '-'])
        #gaps += sum([1 for x in self.bottom if x == '-'])
        if gaps:
            # note: should really implement this - getting an error from a
            # print line is stupid
            msg1 = 'No string coercion method for gapped / single stranded '
            msg2 = 'sequences.'
            print msg1 + msg2
            return NotImplemented
        else:
            return self.top

    def __len__(self):
        '''
        Return length of all DNA (including gaps) in object when built-in
        len function is used.
        '''
        return len(self.top)

    def __add__(self, other):
        '''
        Defines adding with + for DNA objects. Only works for ungapped dsDNA
        and top-only ssDNA.
        '''

        # Check self for terminal gap at the right:
        self_gaps = (self.top[-1] == '-', self.bottom[1] == '-')
        other_gaps = (other.top[1] == '-', other.bottom[-1] == '-')

        # Only time this should fail is when there's a discontinuity. A
        # discontinuity appears when the first and second entries for the gaps
        # tuples are exactly the opposite of each other
        if self_gaps[0] != other_gaps[0] and self_gaps[1] != other_gaps[1]:
            # TODO: pretty sure this is the wrong way to throw an exception
            # here
            print 'Discontinuity at ends of DNA objects - can\'t add.'
            return NotImplemented

        tops = self.top + other.top
        bottoms = other.bottom + self.bottom

        new_instance = self.copy()
        new_instance.top = tops
        new_instance.bottom = bottoms

        return new_instance


class RestrictionSite(object):
    def __init__(self, recognition_site, cut_site, name=None):
        # TODO: input checking / equivalent
        # FIXME: cut_site is 0-indexed for first entry, 1-index for second
        self.recognition_site = recognition_site  # require DNA object
        # cutsite is indexed to leftmost base of restriction site
        self.cut_site = cut_site  # tuple of where top/bottom strands are cut
        # optional name
        self.name = name

    def __repr__(self):
        site = self.recognition_site
        cut = self.cut_site
        cut_symbols = ('|', '|')
        if cut[0] in range(0, len(site)) and cut[1] in range(0, len(site)):
            top_left = str(site[0:self.cut_site[0]])
            top_right = str(site[self.cut_site[0]:])
            top_w_cut = top_left + cut_symbols[0] + top_right

            bottom_left = site[0:self.cut_site[1]].reverse_complement()
            bottom_left = str(bottom_left)[::-1]
            bottom_right = site[self.cut_site[1]:].reverse_complement()
            bottom_right = str(bottom_right)[::-1]
            bottom_w_cut = bottom_left + cut_symbols[1] + bottom_right
        else:
            # TODO: handle enzymes that cut outside of recognition site
            pass

        return '\n'.join([top_w_cut, bottom_w_cut])
