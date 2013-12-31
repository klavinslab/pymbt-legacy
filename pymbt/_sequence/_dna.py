'''DNA object classes.'''
import collections
import os
import re
import shutil
import subprocess
import tempfile
import pymbt.reaction
import pymbt.seqio
from . import utils
from ._sequence import BaseSequence
from pymbt.constants.genbank import TO_PYMBT
from pymbt.constants.molecular_bio import COMPLEMENTS


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
        :returns: pymbt.sequence.DNA instance.
        :raises: ValueError if an element of `features` isn't of type
                 pymbt.sequence.Feature.
                 ValueError if top and bottom strands have different lengths.
                 ValueError if top and bottom strands are not complementary.

        '''
        # Convert to uppercase, run alphabet check
        super(DNA, self).__init__(dna, 'dna', run_checks=run_checks)
        # Set topology
        self.topology = topology
        # Set features
        if features:
            if all([isinstance(feature, Feature) for feature in features]):
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
                    msg = "Top and bottom strands are difference lengths."
                    raise ValueError(msg)
                r_bottom = utils.reverse_complement(self._bottom, 'dna')
                mismatches = [1 for t, b in zip(self._sequence, r_bottom) if
                              (t != b and not (t == '-' or b == '-'))]
                if mismatches:
                    raise ValueError("Bottom strand doesn't match top strand.")
        elif stranded == 'ss':
            self._bottom = ''.join(['-' for x in self._sequence])
        elif stranded == 'ds':
            self._bottom = utils.reverse_complement(self._sequence, 'dna')
        # Set id
        self.id = id
        # Set name
        self.name = name

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
        features shorter than 5 bp.

        :param other: Another sequence.
        :type other: pymbt.sequence.DNA
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

    def ape(self, ape_path=None):
        '''Open in ApE.'''
        cmd = "ApE"
        if ape_path is None:
            # Check for ApE in PATH
            ape_executables = []
            for path in os.environ["PATH"].split(os.pathsep):
                exepath = os.path.join(path, cmd)
                ape_executables.append(os.access(exepath, os.X_OK))
            if not any(ape_executables):
                raise Exception("Ape not in PATH. Use ape_path kwarg.")
        else:
            cmd = ape_path
        # Check whether ApE exists in PATH
        tmp = tempfile.mkdtemp()
        filename = tmp + '/tmp.ape'
        pymbt.seqio.write_dna(self, filename)
        process = subprocess.Popen([cmd, filename])
        # Block until window is closed
        try:
            process.wait()
            shutil.rmtree(tmp)
        except KeyboardInterrupt:
            shutil.rmtree(tmp)

    def bottom(self):
        """Return the raw string of the Crick (bottom) strand.

        :returns: The Crick strand.
        :rtype: str

        """
        return self._bottom

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely-editable copy of the current sequence.
        :rtype: pymbt.sequence.DNA

        '''
        # Significant performance improvements by skipping alphabet check
        features_copy = [feature.copy() for feature in self.features]
        return type(self)(self._sequence, bottom=self._bottom,
                          topology=self.topology, stranded=self.stranded,
                          features=features_copy, id=self.id, name=self.name,
                          run_checks=False)

    def complement(self):
        '''Complement the bases of the sequence.

        :returns: A base-complemented instance of the current sequence.
        :rtype: pymbt.sequence.DNA

        '''
        copy = self.copy()
        code = dict(COMPLEMENTS['dna'])
        copy._sequence = ''.join([code[base] for base in copy._sequence])
        copy._bottom = ''.join([code[base] for base in copy._bottom])
        # Remove features - they make no sense in complement
        copy.features = []
        return copy

    def circularize(self):
        '''Circularize linear DNA.

        :returns: A circularized version of the current sequence.
        :rtype: pymbt.sequence.DNA

        '''
        # FIXME: this should fail for some cases of overhangs.
        copy = self.copy()
        copy.topology = 'circular'
        return copy

    def endswith(self, seq):
        """Report whether parent sequence ends with a query sequence.

        :param seq: Query sequence.
        :type seq: str or pymbt.sequence.DNA
        :returns: Boolean of whether the top strand ends with the query.
        :rtype: bool

        """
        if self._sequence.endswith(str(seq)):
            return True
        else:
            return False

    def extract(self, name, pure=False):
        '''Extract a feature from the DNA sequence.

        :param name: Name of the feature. Must be unique.
        :type name: str
        :param pure: Turn any gaps in the feature into Ns and remove all other
                     features. If False, just extracts start:stop slice.
        :type pure: bool
        :returns: A subsequence from start to stop of the feature.
        :rtype: pymbt.sequence.DNA
        :raises: ValueError if no feature has `name` or more than one match
                 `name`.

        '''
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
                # Turn gaps into Ns
                for gap in extracted.features[0].gaps:
                    for i in range(*gap):
                        extracted[i] = "n"
            return extracted

    def flip(self):
        '''Flip the DNA - swap the top and bottom strands.

        :returns: Flipped DNA (bottom strand is now top strand, etc.).
        :rtype: pymbt.sequence.DNA

        '''
        copy = self.copy()
        copy._sequence, copy._bottom = copy._bottom, copy._sequence
        return copy

    def is_palindrome(self):
        """Report whether sequence is palindromic.

        :returns: Boolean stating whether sequence is a palindrome.
        :rtype: bool

        """
        return utils.palindrome(self)

    def linearize(self, index=0):
        '''Linearize circular DNA at an index.

        :param index: index at which to linearize.
        :type index: int
        :returns: A linearized version of the current sequence.
        :rtype: pymbt.sequence.DNA
        :raises: ValueError if the input is linear DNA.

        '''
        if self.topology == 'linear':
            raise ValueError('Cannot relinearize linear DNA.')
        copy = self.copy()
        copy = copy[index:] + copy[:index]
        copy.topology = 'linear'
        return copy

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str
        :returns: A list of top and bottom strand indices of matches.
        :rtype: list of lists of indices (ints)
        :raises: ValueError if the pattern is longer than either the input
                 sequence (for linear DNA) or twice as long as the input
                 sequence (for circular DNA).

        '''
        # TODO: If linear, should use the methods in BaseSequence
        if self.topology == 'circular':
            if len(pattern) > 2 * len(self):
                raise ValueError('Pattern too long.')
        else:
            if len(pattern) > len(self):
                raise ValueError('Pattern too long.')

        pattern = str(pattern).upper()
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

        return [top_starts, bottom_starts]

    def mw(self):
        """Calculate the molecular weight.

        :returns: The molecular weight of the current sequence.
        :rtype: float

        """
        counter = collections.Counter(self._sequence + self._bottom)
        mw_a = counter["a"] * 313.2
        mw_t = counter["t"] * 304.2
        mw_g = counter["g"] * 289.2
        mw_c = counter["c"] * 329.2
        return mw_a + mw_t + mw_g + mw_c

    def orient(self, index):
        """Orient DNA to index (only applies to circular DNA).

        :param index: DNA position at which to re-zero the DNA.
        :type index: int
        :returns: The current sequence reoriented at `index`.
        :rtype: pymbt.sequence.DNA
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        """
        if self.topology == "linear" and index != 0:
            raise ValueError("Can't reorient linear DNA")
        if index < 0:
            raise ValueError("Reorientation index must be positive")
        else:
            return self[index:] + self[0:index]

    def orient_by_feature(self, featurename):
        """Reorient the DNA based on a feature it contains (circular DNA only).

        :param featurename: A uniquely-named feature.
        :type featurename: str
        :returns: The current sequence reoriented at the start index of a
                  unique feature matching `featurename`.
        :rtype: pymbt.sequence.DNA
        :raises: ValueError if there is no feature of `featurename` or
                 more than one feature matches `featurename`.

        """
        # REFACTOR: Parts are redundant with .extract()
        matched = []
        for feature in self.features:
            if feature.name == featurename:
                matched.append(feature.copy())
        count = len(matched)
        if count == 1:
            return self.orient(matched[0].start)
        elif count > 1:
            raise ValueError("More than one feature has that name.")
        else:
            raise ValueError("No such feature in the sequence.")

    def reverse(self):
        '''Reverse the sequence.

        :returns: A reversed instance of the current sequence.
        :rtype: pymbt.sequence.DNA

        '''
        copy = self.copy()
        copy._sequence = self._sequence[::-1]
        copy._bottom = self._bottom[::-1]
        # Remove features - reversed ones make no sense
        copy.features = []
        return copy

    def reverse_complement(self):
        '''Reverse complement the DNA.

        :returns: A reverse-complemented instance of the current sequence.
        :rtype: pymbt.sequence.DNA

        '''
        copy = self.copy()
        # Store features - they get removed on reverse/complement
        feature_copy = copy.features
        # Note: if sequence is double-stranded, swapping strand is basically
        # (but not entirely) the same thing - gaps affect accuracy.
        copy = copy.reverse()
        copy = copy.complement()
        copy.features = feature_copy

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

    def set_stranded(self, stranded):
        '''Change DNA strandedness

        :param stranded: 'ss' or 'ds' (DNA).
        :type stranded: str
        :returns: The current sequence, converted to ssDNA or dsDNA.
        :rtype: pymbt.sequence.DNA
        :raises: ValueError if `stranded` is not \"ss\" or \"ds\".

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
            if all([char == '-' for char in self._sequence]):
                copy._sequence = reverse_seq._bottom
            elif all([char == '-' for char in self._bottom]):
                copy._bottom = reverse_seq._sequence
            copy.stranded = 'ds'
        else:
            raise ValueError("'stranded' must be 'ss' or 'ds'.")

        return copy

    def startswith(self, seq):
        """Report whether parent sequence starts with a query sequence.

        :param seq: Query sequence.
        :type seq: str or pymbt.sequence.DNA
        :returns: Boolean of whether the top strand starts with the query.
        :rtype: bool

        """
        if self._sequence.startswith(str(seq)):
            return True
        else:
            return False

    def top(self):
        """Return the raw string of the Watson (top) strand.

        :returns: The Watson strand.
        :rtype: str

        """
        return self._sequence

    def transcribe(self):
        '''Transcribe into RNA.

        :returns: An RNA sequence transcribed from the current DNA sequence.
        :rtype: pymbt.sequence.RNA

        '''
        return pymbt.reaction.transcribe(self)

    def _features_on_slice(self, key):
        '''Process features when given a slice (__getitem__). Features that
        would be truncated are simply removed.

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

    def _remove_end_gaps(self):
        """Removes double-stranded gaps from ends of the sequence.

        :returns: The current sequence wiht terminal double-strand gaps ('-')
                  removed.
        :rtype: pymbt.sequence.DNA

        """
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

    def __add__(self, other):
        '''Add DNA together.

        :param other: instance to be added to.
        :type other: compatible sequence object (currently only DNA).
        :returns: Concatenated DNA sequence.
        :rtype: pymbt.sequence.DNA
        :raises: Exception if either sequence is circular.
                 Exception if concatenating a sequence with overhangs would
                 create a discontinuity.

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

    def __contains__(self, query):
        """Defines `query in sequence` operator.

        :param query: query string or DNA sequence
        :type query: str or pymbt.sequence.DNA

        """
        return super(DNA, self).__contains__(query, "n")

    def __delitem__(self, index):
        '''Delete sequence at an index.

        :param index: index to delete
        :type index: int
        :returns: The current sequence with the base at `index` removed.
        :rtype: pymbt.sequence.DNA

        '''
        if self.features:
            self.features = [feature for feature in self.features if index not
                             in range(feature.start, feature.stop)]
            for feature in self.features:
                if feature.start >= index:
                    feature.move(-1)
        super(DNA, self).__delitem__(index)
        bottom_list = list(self._bottom[::-1])
        del bottom_list[index]
        self._bottom = ''.join(bottom_list)[::-1]

    def __getitem__(self, key):
        '''Index and slice sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object
        :returns: A subsequence matching the slice (`key`).
        :rtype: pymbt.sequence.DNA

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

    def __eq__(self, other):
        '''Define equality - sequences, topology, and strandedness are the
        same.

        :returns: Whether current sequence's (Watson and Crick), topology,
                  and strandedness are equivalent to those of another sequence.
        :rtype: bool

        '''
        tops_equal = self._sequence == other._sequence
        bottoms_equal = self._bottom == other._bottom
        topology_equal = self.topology == other.topology
        stranded_equal = self.stranded == other.stranded
        if tops_equal and bottoms_equal and topology_equal and stranded_equal:
            return True
        else:
            return False

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

    def __setitem__(self, index, new_value):
        '''Sets value at index to new value.

        :param index: The index at which the sequence will be modified.
        :type index: int
        :param new_value: The new value at that index
        :type new_value: str or pymbt.sequence.DNA
        :returns: The current sequence with the sequence at `index` replaced
                  with `new_value`.
        :rtype: pymbt.sequence.DNA
        :raises: ValueError if `new_value` is '-'.

        '''
        new_value = str(new_value)
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
        :returns: instance of pymbt.sequence.RestrictionSite

        '''
        self.recognition_site = recognition_site  # require DNA object
        # cutsite is indexed to leftmost base of restriction site
        self.cut_site = cut_site  # tuple of where top/bottom strands are cut
        # optional name
        self.name = name

    def is_palindrome(self):
        '''Report whether sequence is palindromic.

        :returns: Whether the restriction site is a palindrome.
        :rtype: bool

        '''
        return self.recognition_site.is_palindrome()

    def cuts_outside(self):
        '''Report whether the enzyme cuts outside its recognition site.
        Cutting at the very end of the site returns True.

        :returns: Whether the enzyme will cut outside its recognition site.
        :rtype: bool

        '''
        for index in self.cut_site:
            if index < 0 or index > len(self.recognition_site) + 1:
                return True
        return False

    def copy(self):
        '''Return copy of the restriction site.

        :returns: A safely editable copy of the current restriction site.
        :rtype: pymbt.sequence.RestrictionSite

        '''
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
        '''Defines len operator.

        :returns: Length of the recognition site.
        :rtype: int

        '''
        return len(self.recognition_site)


class Primer(object):
    '''A DNA primer - ssDNA with tm, anneal, and optional overhang.'''
    def __init__(self, anneal, tm, overhang=None, name="", note=""):
        '''
        :param anneal: Annealing sequence
        :type anneal: pymbt.sequence.DNA
        :param overhang: Overhang sequence
        :type overhang: pymbt.sequence.DNA
        :param tm: melting temperature
        :type tm: float
        :param name: Optional name of the primer. Used when writing to csv with
                     seqio.write_primers.
        :type name: str
        :param note: Optional description to associate with the primer. Used
                     when writing to csv with seqio.write_primers.
        :type note: str
        :returns: pymbt.sequence.Primer instance.

        '''
        self.tm = tm
        self.anneal = anneal.set_stranded('ss')
        if overhang is not None:
            self.overhang = overhang.set_stranded('ss')
        else:
            self.overhang = DNA('', stranded='ss')
        self.name = name
        self.note = note

    def primer(self):
        '''Produce full (overhang + annealing sequence) primer sequence.

        :returns: The DNA sequence of the primer.
        :rtype: pymbt.sequence.DNA

        '''
        return self.overhang + self.anneal

    def __repr__(self):
        '''Representation of a primer.'''
        if self.overhang:
            return 'Primer: {} Tm: {:.2f}'.format(self.overhang.top().upper() +
                                                  self.anneal.top(), self.tm)
        else:
            return 'Primer: {} Tm: {:.2f}'.format(self.anneal.top(), self.tm)

    def __str__(self):
        '''Coerce DNA object to string.

        :returns: A string of the full primer sequence.
        :rtype: str

        '''
        return str(self.primer())

    def __eq__(self, other):
        '''Define equality - sequences, topology, and strandedness are the
        same.

        :returns: Whether two primers have the same overhang and annealing
                  sequence.
        :rtype: bool

        '''
        anneal_equal = self.anneal == other.anneal
        overhang_equal = self.overhang == other.overhang
        if anneal_equal and overhang_equal:
            return True
        else:
            return False

    def __len__(self):
        '''Define len operator.

        :returns: The length of the full primer sequence.
        :rtype: int

        '''
        return len(self.primer())


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
        :returns: pymbt.sequence.Feature instance.
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
        :rtype: pymbt.sequence.Feature

        '''
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