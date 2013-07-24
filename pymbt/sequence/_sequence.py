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

        '''
        pattern = str(pattern).lower()
        re_pattern = '(?=' + pattern + ')'
        return [index.start() for index in
                re.finditer(re_pattern, self._sequence)]

    def copy(self):
        '''Create a copy of the current instance.'''
        # Significant performance improvements by skipping alphabet check
        return BaseSequence(self._sequence, self._material, run_checks=False)

    def __getitem__(self, key):
        '''Indexing and slicing of sequences.

        :param key: int or slice for subsetting.
        :type key: int or slice

        '''
        copy = self.copy()
        copy._sequence = self._sequence[key]
        return copy

    def __delitem__(self, index):
        '''Deletes sequence at index.

        param index: Index to delete
        type index: int

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
        '''Cast to string.'''
        return self._sequence

    def __len__(self):
        '''Calculate sequence length.'''
        return len(self._sequence)

    def __add__(self, other):
        '''Defines addition.

        :param other: Instance with which to sum.
        :type other: pymbt.sequence._sequence.BaseSequence

        '''
        return BaseSequence(self._sequence + other._sequence, self._material,
                            run_checks=False)

    def __radd__(self, other):
        '''Add unlike types (enables sum function).

        :param other: Object of any other type.
        :param other: pymbt.sequence._sequence.BaseSequence

        '''
        if other == 0 or other is None:
            # For compatibility with sum()
            return self
        elif type(self) != type(other):
            raise TypeError("Can't add {} to {}".format(self, other))
        return self + other

    def __mul__(self, multiplier):
        '''Concatenate copies of the sequence.

        :param multiplier: Factor by which to multiply the sequence.
        :type multiplier: int

        '''
        # Input checking
        if multiplier != int(multiplier):
            raise TypeError("Multiplication by non-integer.")
        return sum(x for x in _decompose(self, multiplier))

    def __eq__(self, other):
        '''Define == operator.

        :param other: Other sequence.
        :type other: pymbt.sequence._sequence.BaseSequence

        '''
        if vars(self) == vars(other):
            return True
        else:
            return False

    def __ne__(self, other):
        '''Define != operator.

        :param other: Other sequence.
        :type other: pymbt.sequence._sequence.BaseSequence

        '''
        try:
            return not (self == other)
        except TypeError:
            return False

    def __contains__(self, pattern):
        '''`x in y`.

        :param pattern: Pattern to find.
        :type pattern: str

        '''
        if str(pattern) in str(self):
            return True
        else:
            return False


def _decompose(string, n):
    '''Given string and multiplier, find 2**n decomposition.

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


class RNA(BaseSequence):
    '''RNA sequence.'''
    def __init__(self, rna, run_checks=True):
        super(RNA, self).__init__(rna, 'rna', run_checks=run_checks)

    def reverse_complement(self):
        '''Reverse complement sequence.'''
        new_instance = self.copy()
        new_instance._sequence = utils.reverse_complement(self._sequence,
                                                          'rna')
        return new_instance


class Peptide(BaseSequence):
    '''Peptide sequence.'''
    def __init__(self, peptide, run_checks=True):
        super(Peptide, self).__init__(peptide, 'peptide',
                                      run_checks=run_checks)
