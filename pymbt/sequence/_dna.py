'''DNA object classes.'''
import re
from pymbt.constants import genbank
from pymbt.sequence import utils
from pymbt.sequence._sequence import BaseSequence
import pymbt.seqio
import tempfile
import shutil
import subprocess


class DNA(BaseSequence):
    '''DNA sequence.'''
    def __init__(self, dna, bottom=None, topology='linear', stranded='ds',
                 features=None, run_checks=True, id=None, name=''):
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
        # Convert to lowercase, run alphabet check
        super(DNA, self).__init__(dna, 'dna', run_checks=run_checks)
        # Set topology
        self.topology = topology
        # Set features
        if features:
            if all(isinstance(feature, Feature) for feature in features):
                self.features = features
            else:
                raise ValueError("non-Feature input for 'features'.")
        else:
            self.features = []
        # Set strandedness
        self.stranded = stranded
        # If bottom was specified, check it + add it
        if bottom:
            self._bottom = bottom
            if run_checks:
                self._bottom = utils.process_seq(bottom, 'dna')
                if len(self._bottom) != len(self._sequence):
                    raise ValueError('Bottom strand is too short.')
                r_bottom = utils.reverse_complement(self._bottom, 'dna')
                mismatches = [1 for t, b in zip(self._sequence, r_bottom) if
                              (t != b and not (t == '-' or b == '-'))]
                if mismatches:
                    raise ValueError("Bottom strand doesn't match top strand.")
        elif stranded == 'ss':
            self._bottom = ''.join('-' for x in self._sequence)
        elif stranded == 'ds':
            self._bottom = utils.reverse_complement(self._sequence, 'dna')
        # Set id
        self.id = id
        # Set name
        self.name = name

    def copy(self):
        '''Create a copy of the current instance.'''
        # Significant performance improvements by skipping alphabet check
        features_copy = [feature.copy() for feature in self.features]
        return type(self)(self._sequence, bottom=self._bottom,
                          topology=self.topology, stranded=self.stranded,
                          features=features_copy, id=self.id, name=self.name,
                          run_checks=False)

    def reverse_complement(self):
        '''Reverse complement the DNA.'''
        copy = self.copy()
        if self.stranded == 'ds':
            copy._sequence = self._bottom
            copy._bottom = self._sequence
        else:
            copy._sequence = utils.reverse_complement(self._bottom, 'dna')
            copy._bottom = utils.reverse_complement(self._sequence, 'dna')

        # Fix features (invert)
        for feature in copy.features:
            # Swap strand
            if feature.strand == 1:
                feature.strand = 0
            else:
                feature.strand = 1
            # Swap start and stop
            feature.start, feature.stop = (feature.stop, feature.start)
            # Adjust start/stop to feature len
            feature.start = len(copy) - feature.start
            feature.stop = len(copy) - feature.stop

        return copy

    def circularize(self):
        '''Circularize linear DNA.'''
        copy = self.copy()
        copy.topology = 'circular'
        return copy

    def linearize(self, index=0):
        '''Linearize circular DNA at an index.

        :param index: index at which to linearize.
        :type index: int

        '''
        if self.topology == 'linear':
            raise ValueError('Cannot relinearize linear DNA.')
        copy = self.copy()
        copy = copy[index:] + copy[:index]
        copy.topology = 'linear'
        return copy

    def set_stranded(self, stranded):
        '''Change DNA strandedness

        :param stranded: 'ss' or 'ds' (DNA).
        :type stranded: str

        '''
        copy = self.copy()
        # Do nothing if already set
        if stranded == self.stranded:
            return copy

        if stranded == 'ss':
            copy._bottom = '-' * len(copy)
            copy.stranded = 'ss'
        elif stranded == 'ds':
            # Find strand that's all gaps (if ss this should be the case)
            reverse_seq = self.reverse_complement()
            if all(char == '-' for char in self._sequence):
                copy._sequence = reverse_seq._sequence
            elif all(char == '-' for char in self._bottom):
                copy._bottom = reverse_seq._bottom
            copy.stranded = 'ds'
        else:
            raise ValueError("'stranded' must be 'ss' or 'ds'.")

        return copy

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str

        '''
        if self.topology == 'circular':
            if len(pattern) > 2 * len(self):
                raise Exception('Pattern longer than 2x (circular) sequence.')
        else:
            if len(pattern) > len(self):
                raise Exception('Pattern longer than (linear) sequence.')

        pattern = str(pattern).lower()
        regex = '(?=' + pattern + ')'

        if self.topology == 'circular':
            r = len(pattern) - 1
            l = len(self) - r + 1
            top = self._sequence[l:] + self._sequence + self._sequence[:r]
            bottom = self._bottom[l:] + self._bottom + self._bottom[:r]
        else:
            top = self._sequence
            bottom = self._bottom

        top_starts = [index.start() for index in re.finditer(regex, top)]
        bottom_starts = [index.start() for index in re.finditer(regex, bottom)]

        # Adjust indices if doing circular search
        if self.topology == 'circular' and len(pattern) > 1:
            top_starts = [start - r + 1 for start in top_starts]
            bottom_starts = [start - r + 1 for start in bottom_starts]

        return (top_starts, bottom_starts)

    def top(self):
        return self._sequence

    def bottom(self):
        return self._bottom

    def extract(self, name):
        '''Extract a feature from the DNA sequence.

        :param name: Name of the feature. Must be unique.
        :type name: str

        '''
        found = [feature for feature in self.features if
                 feature.name == name]
        if not found:
            raise ValueError("Feature list has no feature '{}'".format(name))
        elif len(found) > 1:
            msg = 'Feature name was not unique, found more than one.'
            raise ValueError(msg)
        else:
            return self[found[0].start:found[0].stop]

    def _remove_end_gaps(self):
        '''Removes double-stranded gaps from ends of the sequence.'''
        # TODO: move this to _resect module
        top = self._sequence
        bottom_rev = self._bottom[::-1]

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
        self._sequence = top
        self._bottom = bottom

    def is_palindrome(self):
        '''Report whether sequence is palindromic.'''
        return utils.palindrome(self)

    def ape(self):
        '''Open in ApE.'''
        tmp = tempfile.mkdtemp()
        filename = tmp + '/tmp.ape'
        pymbt.seqio.write_dna(self, filename)
        process = subprocess.Popen(['ApE', filename])
        # Block until window is closed
        try:
            process.wait()
            shutil.rmtree(tmp)
        except KeyboardInterrupt:
            shutil.rmtree(tmp)

    def flip(self):
        '''Flip the DNA - swap the top and bottom strands.'''
        copy = self.copy()
        copy._sequence, copy._bottom = copy._bottom, copy._sequence
        return copy

    def reorient(self, index):
        '''Reorient DNA to index'''
        if index < 0:
            raise ValueError("Reorientation index must be positive")
        else:
            return self[index:] + self[0:index]

    def startswith(self, seq):
        if self._sequence.startswith(str(seq)):
            return True
        else:
            return False

    def endswith(self, seq):
        if self._sequence.endswith(str(seq)):
            return True
        else:
            return False

    def _features_on_slice(self, key):
        '''Process features when given a slice (__getitem__).

        :param key: input to __getitem__, is a slice or int
        :type key: slice or int

        '''
        remaining = []
        if isinstance(key, slice):
            # If a slice, remove stuff that isn't in the slide and adjust
            # feature starts/stops
            if key.step == 1 or key.step is None:
                starts, stops = zip(*[(feature.start, feature.stop) for feature
                                      in self.features])
                for feature in self.features:
                    # if there's key.start only exclude feature.stop too small
                    if key.start and feature.start < key.start:
                        pass
                    # if there's key.stop only exclude feature.stop too big
                    elif key.stop and feature.stop > key.stop:
                        pass
                    else:
                        # if there was a key.start, move feature
                        feature_copy = feature.copy()
                        if key.start:
                            feature_copy.move(-key.start)
                        remaining.append(feature_copy)
        else:
            # Should just be an index. Remove features not in that index.
            # Could type check / cast to int / raise useful exception
            for feature in self.features:
                if feature.start == feature.stop == key:
                    feature.move(key)
                    remaining.append(feature.copy())
        return remaining

    def __getitem__(self, key):
        '''Index and slice sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object

        '''
        copy = super(DNA, self).__getitem__(key)
        copy._bottom = copy._bottom[::-1][key][::-1]
        copy.topology = 'linear'
        if len(copy):
            if copy.features:
                copy.features = self._features_on_slice(key)
        else:
            copy.features = []
        return copy

    def __delitem__(self, index):
        '''Delete sequence at an index.

        param index: index to delete
        type index: int

        '''
        if self.features:
            self.features = [feature for feature in self.features if index not
                             in range(feature.start, feature.stop)]
            for feature in self.features:
                if feature.start >= index:
                    feature.move(-1)
        super(DNA, self).__delitem__(index)
        bottom_list = list(self._bottom[::-1])
        bottom_list.pop(index)
        self._bottom = ''.join(bottom_list)[::-1]

    def __setitem__(self, index, new_value):
        '''Sets value at index to new value.'''
        if new_value == '-':
            raise ValueError("Can't insert gap - split sequence instead.")
        if self.features:
            for i, feature in enumerate(self.features[::-1]):
                if index in range(feature.start, feature.stop):
                    self.features.pop(-(i + 1))
        # setitem on top strand
        super(DNA, self).__setitem__(index, new_value)
        # setitem on bottom strand
        if self.stranded == 'ds':
            sequence_list = list(self._bottom)[::-1]
            sequence_list[index] = str(DNA(new_value).reverse_complement())
            self._bottom = ''.join(sequence_list[::-1])
        else:
            self._bottom = '-' * len(self)

    def __repr__(self):
        '''String to print when object is called directly.'''
        parent = super(DNA, self).__repr__()
        display_bases = 40
        if len(self._sequence) < 90:
            bottom = self._bottom[::-1]
        else:
            rev_bottom = self._bottom[::-1]
            bottom = ''.join([rev_bottom[0:display_bases], ' ... ',
                              rev_bottom[-display_bases:]])
        first_line = '{} {}DNA:'.format(self.topology, self.stranded)
        to_print = '\n'.join([first_line, parent, bottom])
        return to_print

    def __add__(self, other):
        '''Add DNA together.

        :param other: instance to be added to.
        :type other: compatible sequence object (currently only DNA).

        '''
        if self.topology == 'circular' or other.topology == 'circular':
            raise Exception('Can only add linear DNA.')

        discontinuity = [False, False]
        if len(self) != 0 and len(other) != 0:
        # If either is empty, let things proceed anyways
            discontinuity[0] = (self._sequence[-1] == '-' and
                                other._bottom[-1] == '-')
            discontinuity[1] = (self._bottom[0] == '-' and
                                other._sequence[0] == '-')

        for_discontinuity = discontinuity[0]
        rev_discontinuity = discontinuity[1]

        if for_discontinuity or rev_discontinuity:
            msg = "Concatenated DNA would be discontinuous."
            raise Exception(msg)

        if self.stranded == 'ds' or other.stranded == 'ds':
            stranded = 'ds'
        else:
            stranded = 'ss'

        tops = self._sequence + other._sequence
        bottoms = other._bottom + self._bottom
        self_features = [feature.copy() for feature in self.features]
        other_features = [feature.copy() for feature in other.features]
        for feature in other_features:
            feature.move(len(self))

        new_instance = DNA(tops, bottom=bottoms, topology='linear',
                           stranded=stranded, run_checks=False,
                           features=self_features + other_features)

        return new_instance

    def __eq__(self, other):
        '''Define equality - sequences, topology, and strandedness are the
        same.'''
        tops_equal = self._sequence == other._sequence
        bottoms_equal = self._bottom == other._bottom
        topology_equal = self.topology == other.topology
        stranded_equal = self.stranded == other.stranded
        if tops_equal and bottoms_equal and topology_equal and stranded_equal:
            return True
        else:
            return False


class RestrictionSite(object):
    '''Recognition site and properties of a restriction endonuclease.'''
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
        self.recognition_site = recognition_site  # require DNA object
        # cutsite is indexed to leftmost base of restriction site
        self.cut_site = cut_site  # tuple of where top/bottom strands are cut
        # optional name
        self.name = name

    def is_palindrome(self):
        '''Report whether sequence is palindromic.'''
        return self.recognition_site.is_palindrome()

    def cuts_outside(self):
        '''Report whether the enzyme cuts outside its recognition site.

        Cutting at the very end of the site returns True

        '''
        for index in self.cut_site:
            if index < 0 or index > len(self.recognition_site) + 1:
                return True
        return False

    def copy(self):
        '''Return copy of the restriction site.'''
        return RestrictionSite(self.recognition_site, self.cut_site,
                               self.name)

    def __repr__(self):
        '''Represent a restriction site.'''
        site = self.recognition_site
        cut_symbols = ('|', '|')
        if not self.cuts_outside():
            top_left = str(site[0:self.cut_site[0]])
            top_right = str(site[self.cut_site[0]:])
            top_w_cut = top_left + cut_symbols[0] + top_right

            bottom_left = site[0:self.cut_site[1]].reverse_complement()
            bottom_left = str(bottom_left)[::-1]
            bottom_right = site[self.cut_site[1]:].reverse_complement()
            bottom_right = str(bottom_right)[::-1]
            bottom_w_cut = bottom_left + cut_symbols[1] + bottom_right
        else:
            return '\n'.join([site.top() + ' {}'.format(self.cut_site),
                              site.bottom()])

        return '\n'.join([top_w_cut, bottom_w_cut])

    def __len__(self):
        '''Defines len operator.'''
        return len(self.recognition_site)


class Primer(object):
    '''A DNA primer - ssDNA with tm, anneal, and optional overhang.'''
    def __init__(self, anneal, tm, overhang=None):
        '''
        :param anneal: Annealing sequence
        :type anneal: pymbt.sequence.DNA
        :param overhang: Overhang sequence
        :type overhang: pymbt.sequence.DNA
        :param tm: melting temperature
        :type tm: float

        '''
        self.tm = tm
        self.anneal = anneal.set_stranded('ss')
        if overhang is not None:
            self.overhang = overhang.set_stranded('ss')
        else:
            self.overhang = DNA('', stranded='ss')

    def primer(self):
        '''Retrieve full primer sequence.'''
        return self.overhang + self.anneal

    def __repr__(self):
        '''Representation of a primer.'''
        if self.overhang:
            return 'Primer: {} Tm: {:.2f}'.format(self.overhang.top().upper() +
                                                  self.anneal.top(), self.tm)
        else:
            return 'Primer: {} Tm: {:.2f}'.format(self.anneal.top(), self.tm)

    def __str__(self):
        '''Coerce DNA object to string.'''
        return str(self.primer())

    def __eq__(self, other):
        '''Define equality - sequences, topology, and strandedness are the
        same.'''
        anneal_equal = self.anneal == other.anneal
        overhang_equal = self.overhang == other.overhang
        if anneal_equal and overhang_equal:
            return True
        else:
            return False

    def __len__(self):
        '''Define len operator.'''
        return len(self.primer())


class Feature(object):
    '''Represent A DNA feature - annotate and extract sequence by metadata.'''
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

        allowed_types = genbank.TO_PYMBT.keys()

        if feature_type in allowed_types:
            self.feature_type = feature_type
        else:
            msg1 = 'feature_type'
            msg2 = 'must be one of the following: {}'.format(allowed_types)
            raise ValueError(msg1 + msg2)

    def move(self, bases):
        '''Move the start and stop positions.

        :param bases: bases to move - can be negative
        :type bases: int

        '''
        self.start += bases
        self.stop += bases

    def copy(self):
        '''Return a copy of the Feature.'''
        return Feature(self.name, self.start, self.stop, self.feature_type,
                       self.strand)

    def __repr__(self):
        '''Represent a feature.'''
        if self.modified:
            part1 = "(Modified) {} '{}' feature ".format(self.name,
                                                         self.feature_type)
        else:
            part1 = "{} '{}' feature ".format(self.name, self.feature_type)
        part2 = '({0} to {1}) on strand {2}'.format(self.start, self.stop,
                                                    self.strand)
        return part1 + part2

    def __eq__(self, other):
        '''Define equality.'''
        # name is the same
        name_equal = self.name == other.name
        # feature_type is the same
        feature_type_equal = self.feature_type == other.feature_type
        # FIXME: length of feature can't be deduced over origin of circular
        # sequence. Features are very thin so this may not matter.
        # Features are NOT parts and have no associated sequences yet.
        if name_equal and feature_type_equal:
            return True
        else:
            return False

    def __ne__(self, other):
        '''Define inequality.'''
        if not self == other:
            return True
        else:
            return False


def _decompose(string, n):
    '''Given string and multiplier, find n**2 decomposition.

    :param string: input string
    :type string: str
    :param n: multiplier
    :type n: int

    '''
    binary = [int(x) for x in bin(n)[2:]]
    new_string = string
    counter = 1
    while counter <= len(binary):
        if binary[-counter]:
            yield new_string
        new_string += new_string
        counter += 1
