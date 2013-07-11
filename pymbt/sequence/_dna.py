'''
DNA object classes.

'''

import re
from pymbt.sequence import utils

# TODO: figure out what to do with bottom-strand-only ssDNA. Should it be
# flipped automatically or represented as-is? __setitem__ depends on this and
# is currently incomplete. Also has implications for __add__.
#   It should not be flipped automatically (unexpected, limits expressibility)
# TODO: Feature mode is weak, may need to be redesigned
#       Examples:
#           How to track modifications to features? If I have gfp and linearize
#           in it, then recircularize, how can it be restored to 'unmodified'
#           version? Would need sequence information / a history of
#           modifications.
#       Feature updater - if a feature has all negative indices, should update
# TODO: set method for topology


class DNA(object):
    '''
    Core DNA sequence object.

    '''

    def __init__(self, seq, bottom=None, topology='linear', stranded='ds',
                 features=None, run_checks=True, id=None, name=None):
        '''
        :param seq: Input sequence (DNA).
        :type seq: str
        :param bottom: Manual input of bottom-strand sequence. Enables both
                       mismatches and initializing ssDNA.
        :type bottom: str
        :param topology: Topology of DNA - 'linear' or 'circular'.
        :type topology: str
        :param stranded: Strandedness of DNA - 'ss' for single-stranded or
                         'ds' for double-stranded.
        :type stranded: str
        :param features: List of annotated features.
        :type features: list
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :param id: An optional (unique) id field for your DNA sequence.
        :type id: str
        :param name: Optional name field for your DNA sequence.
        :type name: str

        '''

        if run_checks:
            self.top = utils.process_seq(seq, 'dna')
        else:
            self.top = seq

        self.topology = topology
        if not features:
            self.features = []
        else:
            if isinstance(features, Feature):
                self.features = [features]
            elif type(features) == list:
                if all(isinstance(feature, Feature) for feature in features):
                    self.features = features
                else:
                    raise Exception('Invalid feature input')
            else:
                raise Exception('Invalid feature input')
        # TODO: eliminate this attribute - just let it be implicit / checkable
        # with a method e.g. DNA.stranded()
        self.stranded = stranded

        if bottom:
            self.bottom = bottom
            if run_checks:
                self.bottom = utils.process_seq(bottom, 'dna')
                if len(self.bottom) != len(self.top):
                    msg = 'Manually-specified bottom strand is too short.'
                    raise ValueError(msg)
                r_bottom = utils.reverse_complement(self.bottom, 'dna')
                mismatches = [1 for t, b in zip(self.top, r_bottom) if
                              (t != b and (t != '-' or b != '-'))]
                if mismatches:
                    msg = 'Manually-specified bottom strand has mismatches.'
                    raise ValueError(msg)
        elif stranded == 'ss':
            self.bottom = ''.join('-' for x in self.top)
        elif stranded == 'ds':
            self.bottom = utils.reverse_complement(self.top, 'dna')

        # TODO: What should happen when a DNA sequence with an id is modified?
        self.id = id
        self.name = name

    def reverse_complement(self):
        '''
        Reverse complement top and bottom strands.

        '''

        new_instance = self.copy()

        if self.stranded == 'ds':
            new_instance.top = self.bottom
            new_instance.bottom = self.top
        elif self.stranded == 'ss':
            new_instance.top = utils.reverse_complement(self.top, 'dna')

        return new_instance

    def circularize(self):
        '''
        Circularize linear DNA.

        '''

        new_instance = self.copy()
        new_instance.topology = 'circular'

        return new_instance

    def linearize(self, index=0):
        '''
        Linearize circular DNA at a specific index.

        :param index: index at which to linearize.
        :type index: int

        '''

        if self.topology == 'linear':
            raise Exception('Cannot relinearize linear DNA.')

        new_instance = self.copy()
        new_instance = new_instance[index:] + new_instance[:index]
        new_instance.topology = 'linear'

        return new_instance

    def five_resect(self, n_bases):
        '''
        Remove bases from 5' end of top strand.


        :param n_bases: Number of bases cut back.
        :type n_bases: int

        '''

        new_instance = self.copy()

        new_top = '-' * min(len(self.top), n_bases) + self.top[n_bases:]
        new_instance.top = new_top
        new_instance.remove_end_gaps()
        if n_bases >= len(self):
            self = self.set_stranded('ss')

        return new_instance

    def three_resect(self, n_bases):
        '''
        Remove bases from 3' end of top strand.

        :param n_bases: Number of bases cut back.
        :type n_bases: int

        '''

        # TODO: if you find double end gaps, should make topology linear
        new_instance = self.copy()

        new_top = self.top[:-n_bases] + '-' * min(len(self.top), n_bases)
        new_instance.top = new_top
        new_instance.remove_end_gaps()
        if n_bases >= len(self):
            self = self.set_stranded('ss')

        return new_instance

    def set_stranded(self, stranded):
        '''
        Change DNA strandedness

        :param stranded: 'ss' or 'ds' (DNA).
        :type stranded: str

        '''

        new_instance = self.copy()
        # Do nothing if already set
        if stranded == self.stranded:
            return new_instance

        if stranded == 'ss':
            new_instance.bottom = '-' * len(new_instance)
            new_instance.stranded = 'ss'
        elif stranded == 'ds':
            # Find strand that's all gaps (if ss this should be the case)
            reverse_seq = self.reverse_complement()
            if all(char == '-' for char in self.top):
                new_instance.top = reverse_seq.bottom
            elif all(char == '-' for char in self.bottom):
                new_instance.bottom = reverse_seq.top
            new_instance.stranded = 'ds'
        else:
            raise ValueError("'stranded' must be 'ss' or 'ds'.")

        return new_instance

    def locate(self, pattern):
        '''
        Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str

        '''

        if len(pattern) > 2 * len(self):
            raise Exception('Pattern must be less than 2 x sequence length.')
        if not pattern:
            return [[], []]

        pattern = pattern.lower()
        regex = '(?=' + pattern + ')'

        if self.topology == 'circular':
            roff = len(pattern) - 1
            loff = len(self) - roff + 1
            top = self.top[loff:] + self.top + self.top[0:roff]
            bottom = self.bottom[loff:] + self.bottom + self.bottom[0:roff]
        else:
            top = self.top
            bottom = self.bottom

        top_starts = [index.start() for index in re.finditer(regex, top)]
        bottom_starts = [index.start() for index in re.finditer(regex, bottom)]

        # Adjust indices if doing circular search
        if self.topology == 'circular':
            top_starts = [start - roff + 1 for start in top_starts]
            bottom_starts = [start - roff + 1 for start in bottom_starts]

        return (top_starts, bottom_starts)

    def extract(self, feature_name):
        '''
        Extract a feature from the DNA sequence.

        :param feature_name: Name of the future. Must be unique.
        :type feature_name: str

        '''

        found_feature = False
        for feature in self.features:
            if feature.name == feature_name:
                if found_feature:
                    msg = 'Feature name was not unique, found more than one.'
                    raise ValueError(msg)
                else:
                    found_feature = feature
        if found_feature:
            return self[found_feature.start:found_feature.stop]
        else:
            raise ValueError('Feature name does not appear in features list')

    def copy(self):
        '''
        Create a copy of the current instance.

        '''

        # Alphabet checking disabled on copy to improve performance
        new_instance = DNA(self.top, bottom=self.bottom,
                           topology=self.topology, stranded=self.stranded,
                           features=self.features, run_checks=False,
                           id=self.id, name=self.name)

        return new_instance

    def remove_end_gaps(self):
        '''
        Removes double-stranded gaps from ends of the sequence.

        '''

        top = self.top
        bottom_rev = self.bottom[::-1]

        top_lstrip = top.lstrip('-')
        bottom_lstrip = bottom_rev.lstrip('-')
        lstrip_len = max(len(top_lstrip), len(bottom_lstrip))
        top = top[-lstrip_len:]
        bottom_rev = bottom_rev[-lstrip_len:]

        top_rstrip = top.rstrip('-')
        bottom_rstrip = bottom_rev.rstrip('-')
        rstrip_len = max(len(top_rstrip), len(bottom_rstrip))
        top = top[0:rstrip_len]
        bottom_rev = bottom_rev[0:rstrip_len]

        bottom = bottom_rev[::-1]
        self.top = top
        self.bottom = bottom

    def __getitem__(self, key):
        '''
        Indexing and slicing of sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object

        '''

        new_instance = self.copy()
        new_instance.top = new_instance.top[key]
        new_instance.bottom = new_instance.bottom[::-1][key][::-1]

        new_instance.topology = 'linear'

        return new_instance

    def __delitem__(self, index):
        '''
        Deletes sequence at index.

        param index: index to delete
        type index: int

        '''

        top_list = list(self.top)
        bottom_list = list(self.bottom[::-1])

        del top_list[index]
        del bottom_list[index]

        self.top = ''.join(top_list)
        self.bottom = ''.join(bottom_list)[::-1]

    def __setitem__(self, index, new_value):
        '''
        Sets index value to new value.

        '''

        if new_value == '-':
            raise ValueError("Can't insert gap - split sequence instead.")

        insert = DNA(str(new_value))
        new_top = insert.top
        if self.stranded == 'ds':
            new_bottom = insert.bottom
        else:
            new_bottom = '-'

        top_list = list(self.top)
        bottom_list = list(self.bottom[::-1])

        top_list[index] = new_top
        bottom_list[index] = new_bottom

        self.top = ''.join(top_list)
        self.bottom = ''.join(bottom_list)[::-1]

    def __repr__(self):
        '''
        String to print when object is called directly.

        '''

        show = 40
        bottom = self.bottom[::-1]

        if len(self.top) < 90:
            top = self.top
            bottom = self.bottom[::-1]
        else:
            top = ''.join([self.top[0:show], ' ... ', self.top[-show:]])
            bottom = ''.join([bottom[0:show], ' ... ', bottom[-show:]])

        first_line = '{} {}DNA:'.format(self.topology, self.stranded)
        to_print = '\n'.join([first_line, top, bottom])

        return to_print

    def __str__(self):
        '''
        Coerce DNA object to string. Only works for ungapped dsDNA
        and top-strand ssDNA.

        '''

        if '-' in self.top:
            msg = 'No string coercion method sequences with top-strand gaps.'
            raise Exception(msg)
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

        :param other: instance to be added to.
        :type other: compatible sequence object (currently only DNA).

        '''

        if self.topology == 'circular' or other.topology == 'circular':
            raise Exception('Can only add linear DNA.')

        discontinuity = [False, False]
        if len(self) != 0 and len(other) != 0:
        # If either is empty, let things proceed anyways
            discontinuity[0] = self.top[-1] == '-' and other.bottom[-1] == '-'
            discontinuity[1] = self.bottom[0] == '-' and other.top[0] == '-'

        for_discontinuity = discontinuity[0]
        rev_discontinuity = discontinuity[1]

        if for_discontinuity or rev_discontinuity:
            msg = "Concatenated DNA would be discontinuous."
            raise Exception(msg)

        if self.stranded == 'ds' or other.stranded == 'ds':
            stranded = 'ds'
        else:
            stranded = 'ss'

        tops = self.top + other.top
        bottoms = other.bottom + self.bottom

        new_instance = DNA(tops, bottom=bottoms, topology='linear',
                           stranded=stranded, run_checks=False)

        return new_instance

    def __radd__(self, other):
        '''
        Add unlike types (enables sum function).

        :param self: object of self type.
        :type self: DNA
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
        Multiply DNA by an integer to create concatenation.

        :param multiplier: Factor by which to multiply the sequence.
        :type multiplier: int

        '''

        # Input checking
        if multiplier != int(multiplier):
            raise TypeError("can't multiply sequence by non-integer.")
        if self.topology == 'circular':
            raise ValueError("Can't multiply circular DNA")

        # Test concatenation by adding once
        try:
            self + self
        except:
            raise Exception('Failed to add, so cannot multiply.')

        # Isolate top and bottom strands, multiply strings, recreate DNA
        tops = self.top * multiplier
        bottoms = self.bottom * multiplier

        new_instance = DNA(tops, bottom=bottoms, topology=self.topology,
                           stranded=self.stranded, run_checks=False)

        return new_instance

    def __eq__(self, other):
        '''
        Test DNA object equality.

        '''

        if vars(self) == vars(other):
            return True
        else:
            return False

    def __ne__(self, other):
        '''
        Test DNA object equality.

        '''

        return not (self == other)

    def __contains__(self, pattern):
        '''
        Specially-defined `x in y` behavior

        '''

        if pattern in self.top or pattern in self.bottom:
            return True
        else:
            return False


class RestrictionSite(object):
    '''
    Defines the recognition site and properties of a restriction endonuclease.

    '''

    def __init__(self, recognition_site, cut_site, name=None):
        '''
        :param recognition_site: Input sequence.
        :type recognition_site: pymbt.sequence.DNA
        :param cut_site: 0-indexed indices where DNA is nicked (top, then
                         bottom strand). For an n-sized recognition site, there
                         are n + 1 positions at which to cut.
        :type cut_site: 2-tuple.
        :param name: Identifier of this restriction site
        :type name: str

        '''
        # TODO: input checking / equivalent
        self.recognition_site = recognition_site  # require DNA object
        # cutsite is indexed to leftmost base of restriction site
        self.cut_site = cut_site  # tuple of where top/bottom strands are cut
        # optional name
        self.name = name

    def __repr__(self):
        '''
        Restriction site representation / string coercion.

        '''

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


class Feature(object):
    '''
    A DNA feature - annotate and extract sequence by metadata.

    '''

    def __init__(self, name, start, stop, feature_type, strand=0):
        '''
        :param name: Name of the feature. Used during feature extraction.
        :type name: str
        :param start: Where the feature starts
        :type start: int
        :param stop: Where the feature stops
        :type stop: int
        :param feature_type: The type of the feature. Allowed types:
                                'coding', 'primer', 'promoter', 'terminator',
                                'rbs'
        :type name: str
        :param strand: Watson (0) or Crick (1) strand of the feature.
        :type strand: int

        '''

        self.name = name
        self.start = int(start)
        self.stop = int(stop)
        self.modified = False
        self.strand = strand

        allowed_types = ['coding', 'primer', 'promoter', 'terminator', 'rbs',
                         'misc', 'origin']

        if feature_type in allowed_types:
            self.feature_type = feature_type
        else:
            msg1 = 'feature_type'
            msg2 = 'must be one of the following: {}'.format(allowed_types)
            raise ValueError(msg1 + msg2)

    def move(self, bases):
        '''
        Move the start and stop positions.

        :param bases: bases to move - can be negative
        :type bases: int

        '''

        self.start += bases
        self.stop += bases

    def __repr__(self):
        '''
        Representation of the feature.

        '''

        if self.modified:
            part1 = "(Modified) {} '{}' feature ".format(self.name,
                                                         self.feature_type)
        else:
            part1 = "{} '{}' feature ".format(self.name, self.feature_type)
        part2 = '({0} to {1}) on strand {2}'.format(self.start, self.stop,
                                                    self.strand)
        return part1 + part2
