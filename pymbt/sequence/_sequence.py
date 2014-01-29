'''Base sequence classes.'''
import re
from . import utils
from pymbt.constants.genbank import TO_PYMBT
from pymbt.constants.molecular_bio import ALPHABETS, COMPLEMENTS


class BaseSequence(object):
    '''Base sequence class.'''
    def __init__(self, sequence, material, features=None, run_checks=True):
        '''
        :param sequence: Input sequence.
        :type sequence: str
        :param material: Material type (dna, rna, peptide)
        :type material: str
        :param features: List of annotated features.
        :type features: list
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :returns: pymbt.sequence.BaseSequence instance.

        '''
        if run_checks:
            self._sequence = utils.process_seq(sequence, material)
        else:
            self._sequence = sequence
        self._material = material

        # Set features
        if features:
            # FIXME: violated duck-typing here
            if all([isinstance(feature, Feature) for feature in features]):
                self.features = features
            else:
                raise ValueError("non-Feature input for 'features'.")
        else:
            self.features = []

    def annotate_from_library(self, library, wipe=True, shortest=6):
        """Annotate the sequence using the features of another. Ignores
        features shorter than 5 bp.

        :param library: list of features with a .sequence attribute
        :type library: list of features with a .sequence attribute
        :param wipe: Remove (wipe) the current features first.
        :type wipe: bool
        :param shortest: Features shorters than this will be ignored.
        :type shortest: int

        """
        copy = self.copy()
        if wipe:
            copy.features = []

        # Make sure features are unique and above 'shortest' parameter
        unique_features = []
        unique_sequences = []
        for feature in library:
            if len(feature.sequence) >= shortest:
                if feature.sequence not in unique_sequences:
                    unique_features.append(feature)
                    unique_sequences.append(feature.sequence)
        # Match features
        for feature in unique_features:
            if feature.sequence in copy:
                match_location = copy.locate(feature.sequence)
                # Only annotate the top strand if sequence is palindrome
                if feature.sequence.is_palindrome():
                    match_location[1] = []
                # Which strand is it on?
                for i, strand in enumerate(match_location):
                    for match in strand:
                        new_feature = feature.copy()
                        length = new_feature.stop - new_feature.start
                        if i == 0:
                            # Watson strand
                            new_feature.start = match
                            new_feature.stop = new_feature.start + length
                        else:
                            # Crick strand
                            new_feature.stop = len(copy) - match
                            new_feature.start = new_feature.stop - length
                            # If feature found on other strand, update
                            new_feature.strand = abs(new_feature.strand - 1)
                        # Modulo just in case it spans the origin
                        new_feature.start = new_feature.start % len(copy)
                        new_feature.stop = new_feature.stop % len(copy)
                        copy.features.append(new_feature)
        return copy

    def annotate_from_other(self, other, wipe=True, shortest=6):
        """Annotate the sequence using the features of another. Ignores
        features shorter than 5.

        :param other: Another sequence.
        :type other: pymbt.DNA
        :param wipe: Remove (wipe) the current features first.
        :type wipe: bool
        :param shortest: Features shorters than this will be ignored.
        :type shortest: int

        """
        # Generate feature library
        features = [feature.copy() for feature in other.features]
        for feature in features:
            feature.sequence = other[feature.start:feature.stop]

        # Annotate from library
        return self.annotate_from_library(features, wipe=wipe,
                                          shortest=shortest)

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely editable copy of the current sequence.

        '''
        # Significant performance improvements by skipping alphabet check
        features_copy = [feature.copy() for feature in self.features]
        return type(self)(self._sequence, self._material,
                          features=features_copy, run_checks=False)

    def endswith(self, seq):
        """Report whether parent sequence ends with a query sequence.

        :param seq: Query sequence.
        :returns: Boolean of whether the sequence ends with the query.
        :rtype: bool

        """
        if self._sequence.endswith(str(seq)):
            return True
        else:
            return False

    def extract(self, name, any_char, pure=False):
        '''Extract a feature from the sequence.

        :param name: Name of the feature. Must be unique.
        :type name: str
        :param pure: Turn any gaps in the feature into Ns or Xs and remove all
                     other features. If False, just extracts start:stop slice.
        :type pure: bool
        :returns: A subsequence from start to stop of the feature.
        :raises: ValueError if no feature has `name` or more than one match
                 `name`.

        '''
        # TODO: reconsider 'pure' default (could be True instead)
        found = [feature for feature in self.features if
                 feature.name == name]
        if not found:
            raise ValueError("Feature list has no feature '{}'".format(name))
        elif len(found) > 1:
            msg = 'Feature name was not unique, found more than one.'
            raise ValueError(msg)
        else:
            extracted = self[found[0].start:found[0].stop]
            if pure:
                # Keep only the feature specified
                extracted.features = [found[0]]
                # Turn gaps into Ns or Xs
                for gap in extracted.features[0].gaps:
                    for i in range(*gap):
                        extracted[i] = any_char
            return extracted

    def insert(self, sequence, index):
        '''Insert a sequence at index.

        :param sequence: Sequence to insert
        :param index: (internal) index at which to insert the sequence.
        :type index: int

        '''
        range_error = IndexError("Invalid index - must be between 1 and " +
                                 "length - 1.")
        if index == 0:
            raise range_error
        try:
            self[index]
        except IndexError:
            raise range_error

        return self[0:index] + sequence + self[index:]

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str
        :returns: Indices of pattern matches.
        :rtype: list of ints

        '''
        pattern = str(pattern).upper()
        re_pattern = '(?=' + pattern + ')'
        return [index.start() for index in
                re.finditer(re_pattern, self._sequence)]

    def startswith(self, query):
        """Report whether parent sequence starts with a query sequence.

        :param seq: Query sequence.
        :type seq: str or pymbt.DNA
        :returns: Boolean of whether the top strand starts with the query.
        :rtype: bool

        """
        if self._sequence.startswith(str(query)):
            return True
        else:
            return False

    def __add__(self, other):
        '''Defines addition.

        :param other: Instance with which to sum.
        :type other: pymbt.sequence.BaseSequence
        :returns: Concatenated sequence.
        :rtype: pymbt.sequence.BaseSequence

        '''
        self_features = [feature.copy() for feature in self.features]
        other_features = [feature.copy() for feature in other.features]
        for feature in other_features:
            feature.move(len(self))
        features = self_features + other_features

        return BaseSequence(self._sequence + other._sequence, self._material,
                            features=features, run_checks=False)

    def __contains__(self, query, any_char):
        '''`x in y`.

        :param query: Query (i.e. exact pattern) sequence to find.
        :type query: str
        :param any_char: Character to use for any match (.*).
        :type any_char: str
        :returns: Whether the query is found in the current sequence.
        :rtype: bool

        '''
        query_str = str(query).upper()
        query_str = re.sub(any_char, ".", query_str)
        if re.search(query_str, str(self)):
            return True
        else:
            return False

    def __delitem__(self, index):
        '''Deletes sequence at index.

        :param index: Index to delete
        :type index: int
        :returns: The current sequence with the moiety at `index` removed.
        :rtype: pymbt.sequence.BaseSequence

        '''
        if self.features:
            self.features = [feature for feature in self.features if index not
                             in range(feature.start, feature.stop)]
            for feature in self.features:
                if feature.start >= index:
                    feature.move(-1)

        sequence_list = list(self._sequence)
        del sequence_list[index]
        self._sequence = ''.join(sequence_list)

    def __eq__(self, other):
        '''Define == operator. True if sequences are the same.

        :param other: Other sequence.
        :type other: pymbt.sequence._sequence.BaseSequence
        :returns: Whether two sequences have the same base string (sequence).
        :rtype: bool

        '''
        if self._sequence == other._sequence:
            return True
        else:
            return False

    def __getitem__(self, key):
        '''Indexing and slicing of sequences.

        :param key: int or slice for subsetting.
        :type key: int or slice
        :returns: Slice of the current sequence.
        :rtype: pymbt.sequence.BaseSequence

        '''
        copy = self.copy()
        copy._sequence = self._sequence[key]

        def in_slice(feature):
            if key.start and feature.start < key.start:
                return False
            elif key.stop and feature.stop > key.stop:
                return False
            else:
                return True

        if not len(copy):
            copy.features = []
            return copy

        if copy.features:
            if isinstance(key, slice):
                # If a slice, remove stuff that isn't in the slide and
                # adjust feature starts/stops
                if key.step == 1 or key.step is None:
                    copy.features = [feature.copy() for feature in
                                     self.features if in_slice(feature)]
                    if key.start:
                        for feature in copy.features:
                            feature.move(-key.start)
                else:
                    copy.features = []
            else:
                copy.features = [feature.copy() for feature in self.features if
                                 feature.start == feature.stop == key]
                for feature in copy.features:
                    feature.move(key)
        return copy

    def __len__(self):
        '''Calculate sequence length.

        :returns: The length of the sequence.
        :rtype: int

        '''
        return len(self._sequence)

    def __mul__(self, n):
        '''Concatenate copies of the sequence.

        :param n: Factor by which to multiply the sequence.
        :type n: int
        :returns: The current sequence repeated n times.
        :rtype: pymbt.sequence.BaseSequence
        :raises: TypeError if n is not an integer.

        '''
        # Input checking
        if n != int(n):
            raise TypeError("Multiplication by non-integer.")
        return sum([x for x in _decompose(self, n)])

    def __ne__(self, other):
        '''Define != operator.

        :param other: Other sequence.
        :type other: pymbt.sequence._sequence.BaseSequence
        :returns: The opposite of ==.
        :rtype: bool

        '''
        try:
            return not (self == other)
        except TypeError:
            return False

    def __radd__(self, other):
        '''Add unlike types (enables sum function).

        :param other: Object of any other type.
        :param other: pymbt.sequence._sequence.BaseSequence
        :returns: Concatenated sequence.
        :rtype: pymbt.sequence.BaseSequence

        '''
        if other == 0 or other is None:
            # For compatibility with sum()
            return self
        elif type(self) != type(other):
            raise TypeError("Can't add {} to {}".format(self, other))
        return self + other

    def __repr__(self):
        '''String to print when object is called directly.'''
        display_bases = 40
        if len(self._sequence) < 90:
            sequence = self._sequence
        else:
            sequence = ''.join([self._sequence[:display_bases], ' ... ',
                                self._sequence[-display_bases:]])
        return sequence

    def __setitem__(self, index, new_value):
        '''Sets index value to new value.

        :param index: Index to modify.
        :type index: int
        :param new_value: Value to input.
        :type new_value: str or pymbt.sequence._sequence.BaseSequence
        :returns: The current sequence with the moiety at `index` replaced
                  by `new_value`.
        :rtype: pymbt.sequence.BaseSequence

        '''
        if self.features:
            for i, feature in enumerate(self.features[::-1]):
                if index in range(feature.start, feature.stop):
                    self.features.pop(-(i + 1))

        sequence_list = list(self._sequence)
        sequence_list[index] = str(BaseSequence(new_value, self._material))
        self._sequence = ''.join(sequence_list)

    def __str__(self):
        '''Cast to string.

        :returns: A string of the current sequence
        :rtype: str

        '''
        return self._sequence


class NucleotideSequence(BaseSequence):
    '''Nucleotide sequence class.'''
    def is_palindrome(self):
        """Report whether sequence is palindromic.

        :returns: Boolean stating whether sequence is a palindrome.
        :rtype: bool

        """
        return utils.palindrome(self)


def _decompose(string, n):
    '''Given string and multiplier n, find m**2 decomposition.

    :param string: input string
    :type string: str
    :param n: multiplier
    :type n: int
    :returns: generator that produces m**2 * string if m**2 is a factor of n
    :rtype: generator of 0 or 1

    '''
    binary = [int(x) for x in bin(n)[2:]]
    new_string = string
    counter = 1
    while counter <= len(binary):
        if binary[-counter]:
            yield new_string
        new_string += new_string
        counter += 1


class Feature(object):
    '''Represent A DNA feature - annotate and extract sequence by metadata.'''
    def __init__(self, name, start, stop, feature_type, strand=0, gaps=[]):
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
        :param gaps: Gap locations if the feature has gaps.
        :type gaps: list of coordinates (2-tuple/list)
        :returns: pymbt.Feature instance.
        :raises: ValueError if `feature_type` is not in
                 pymbt.constants.genbank.TO_PYMBT.

        '''
        self.name = name
        self.start = int(start)
        self.stop = int(stop)
        self.modified = False
        self.strand = strand
        self.gaps = gaps

        allowed_types = TO_PYMBT.keys()

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
        '''Return a copy of the Feature.

        :returns: A safely editable copy of the current feature.
        :rtype: pymbt.Feature

        '''
        return type(self)(self.name, self.start, self.stop, self.feature_type,
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
        '''Define equality.

        :returns: Whether the name and feature type are the same.
        :rtype: bool

        '''
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


def reverse_complement(sequence, material):
    """Reverse complement a sequence.

    :param sequence: Sequence to reverse complement
    :type sequence: str
    :param material: dna, rna, or peptide.
    :type material: str
    """
    # TODO: put in _sequence module and import
    code = dict(COMPLEMENTS[material])
    # TODO: see if using reversed() here has a speed cost
    # FIXME: reverse_complement is redundant with flip?
    reverse_sequence = sequence[::-1]
    return ''.join([code[base] for base in reverse_sequence])


def check_alphabet(seq, material):
    '''Verify that a given string is valid DNA, RNA, or peptide characters.

    :param seq: DNA, RNA, or peptide sequence.
    :type seq: str
    :param material: Input material - 'dna', 'rna', or 'pepide'.
    :type sequence: str
    :returns: Whether the `seq` is a valid string of `material`.
    :rtype: bool
    :raises: ValueError if `material` isn't \"dna\", \"rna\", or \"peptide\".
             ValueError if `seq` contains invalid characters for its
             material type.

    '''
    errs = {'dna': 'DNA', 'rna': 'RNA', 'peptide': 'peptide'}
    if material == 'dna' or material == 'rna' or material == 'peptide':
        alphabet = ALPHABETS[material]
        err_msg = errs[material]
    else:
        msg = "Input material must be 'dna', 'rna', or 'peptide'."
        raise ValueError(msg)
    # This is a bottleneck when modifying sequence - hence the run_checks
    # optional parameter in sequence objects..
    # First attempt with cython was slower. Could also try pypy.
    if re.search('[^' + alphabet + ']', seq):
        raise ValueError('Encountered a non-%s character' % err_msg)


def process_seq(seq, material):
    '''Validate and process sequence inputs.

    :param seq: input sequence
    :type seq: str
    :param material: DNA, RNA, or peptide
    :type: str
    :returns: Uppercase version of `seq` with the alphabet checked by
              check_alphabet().
    :rtype: str

    '''
    check_alphabet(seq, material)
    seq = seq.upper()
    return seq


def palindrome(seq):
    '''Test whether a sequence is palindrome.

    :param seq: Sequence to analyze (DNA or RNA).
    :type seq: pymbt.DNA or pymbt.RNA
    :returns: Whether a sequence is a palindrome.
    :rtype: bool

    '''
    seq_len = len(seq)
    if seq_len % 2 == 0:
        # Sequence has even number of bases, can test non-overlapping seqs
        wing = seq_len / 2
        l_wing = seq[0: wing]
        r_wing = seq[wing:]
        if l_wing == r_wing.reverse_complement():
            return True
        else:
            return False
    else:
        # Sequence has odd number of bases and cannot be a palindrome
        return False
