'''Base sequence classes.'''
import re
from pymbt.sequence import utils


class BaseSequence(object):
    '''Base sequence class.'''
    def __init__(self, sequence, material, run_checks=True):
        '''
        :param sequence: Input sequence.
        :type sequence: str
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

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str
        :returns: Indices of pattern matches.
        :rtype: list of ints

        '''
        pattern = str(pattern).lower()
        re_pattern = '(?=' + pattern + ')'
        return [index.start() for index in
                re.finditer(re_pattern, self._sequence)]

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely editable copy of the current sequence.
        :rtype: pymbt.sequence.BaseSequence

        '''
        # Significant performance improvements by skipping alphabet check
        return type(self)(self._sequence, self._material, run_checks=False)

    def __getitem__(self, key):
        '''Indexing and slicing of sequences.

        :param key: int or slice for subsetting.
        :type key: int or slice
        :returns: Slice of the current sequence.
        :rtype: pymbt.sequence.BaseSequence

        '''
        copy = self.copy()
        copy._sequence = self._sequence[key]
        return copy

    def __delitem__(self, index):
        '''Deletes sequence at index.

        :param index: Index to delete
        :type index: int
        :returns: The current sequence with the moiety at `index` removed.
        :rtype: pymbt.sequence.BaseSequence

        '''
        sequence_list = list(self._sequence)
        sequence_list.pop(index)
        self._sequence = ''.join(sequence_list)

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
        sequence_list = list(self._sequence)
        sequence_list[index] = str(BaseSequence(new_value, self._material))
        self._sequence = ''.join(sequence_list)

    def __repr__(self):
        '''String to print when object is called directly.'''
        display_bases = 40
        if len(self._sequence) < 90:
            sequence = self._sequence
        else:
            sequence = ''.join([self._sequence[:display_bases], ' ... ',
                                self._sequence[-display_bases:]])
        return sequence

    def __str__(self):
        '''Cast to string.

        :returns: A string of the current sequence
        :rtype: str

        '''
        return self._sequence

    def __len__(self):
        '''Calculate sequence length.

        :returns: The length of the sequence.
        :rtype: int

        '''
        return len(self._sequence)

    def __add__(self, other):
        '''Defines addition.

        :param other: Instance with which to sum.
        :type other: pymbt.sequence.BaseSequence
        :returns: Concatenated sequence.
        :rtype: pymbt.sequence.BaseSequence

        '''
        return BaseSequence(self._sequence + other._sequence, self._material,
                            run_checks=False)

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

    def __contains__(self, query, any_char):
        '''`x in y`.

        :param query: Query (i.e. exact pattern) sequence to find.
        :type query: str
        :param any_char: Character to use for any match (.*).
        :type any_char: str
        :returns: Whether the query is found in the current sequence.
        :rtype: bool

        '''
        query_str = str(query).lower()
        query_str = re.sub(any_char, ".", query_str)
        if re.search(query_str, str(self)):
            return True
        else:
            return False


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
