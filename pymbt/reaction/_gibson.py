'''Gibson reaction simulation.'''
import pymbt.analysis


class AmbiguousGibsonError(ValueError):
    '''Exception to raise when Gibson is ambiguous.'''
    pass


class Gibson(object):
    '''Gibson reaction.'''
    def __init__(self, seq_list):
        '''
        :param seq_list: list of DNA sequences to Gibson
        :type seq_list: list of pymbt.DNA
        :returns: pymbt.reaction.Gibson instance.
        :raises: ValueError if any input sequences are circular DNA.

        '''
        # FIXME: Preserve features in overlap
        # TODO: set a max length?
        # TODO: add 'expected' keyword argument somewhere to automate
        # validation of output
        self._seq_list = []
        # Remove any redundant (identical) sequences
        for seq in seq_list:
            if seq not in self._seq_list:
                self._seq_list.append(seq)
        if any([seq.topology == 'circular' for seq in self._seq_list]):
            raise ValueError('Input sequences must be linear, not circular.')

    def run_circular(self, homology_min=10, tm_min=65):
        '''Attempt to produce circular fragment from input fragments.

        :param homology_min: minimum bp of homology allowed
        :type homology_min: int
        :returns: Gibson-assembled (circular) DNA.
        :rtype: pymbt.DNA

        '''
        return self._run(linear=False, homology=homology_min, tm=tm_min)

    def run_linear(self, homology_min=10, tm_min=65):
        '''Attempt to produce linear fragment from input fragments.

        :param homology_min: minimum bp of homology allowed
        :type homology_min: int
        :returns: Gibson-assembled (linear) DNA.
        :rtype: pymbt.DNA

        '''
        return self._run(linear=True, homology=homology_min, tm=tm_min)

    def _run(self, linear=False, homology=10, tm=65):
        '''Run the Gibson reaction.

        :param linear: Return a linear fragment rather than attempting to
                       generate a circular construct.
        :type linear: bool
        :param homology_min: minimum bp of homology allowed
        :type homology_min: int
        :returns: Gibson-assembled DNA.
        :rtype: pymbt.DNA

        '''
        # Copy input list
        self._working_list = self._seq_list[:]
        # Attempt to fuse fragments together until only one is left
        while len(self._working_list) > 1:
            self._find_fuse_next(homology, tm)
        if not linear:
            # Fuse the final fragment to itself
            self._fuse_last(homology, tm)
        return self._working_list[0]

    def _find_fuse_next(self, homology, tm):
        '''Find the next sequence to fuse, and fuse it (or raise exception).

        :param homology: length of terminal homology in bp
        :type homology: int
        :raises: AmbiguousGibsonError if there is more than one way for the
                 fragment ends to combine.

        '''
        # 1. Analyze all non-first sequences for matches
        pattern = str(self._working_list[0])
        targets = self._working_list[1:]
        report = [homology_report(pattern, target) for target in targets]

        # 2. Throw out matches that are too short or have too low of tm
        for i, target in enumerate(report):
            for j, side in enumerate(target):
                # If too short, throw out match
                report[i][j] = [x for x in side if x > homology]
                # If tm is too low, throw out match
                if j == 0:
                    seqs = [targets[i][:x] for x in side]
                    report[i][j] = [x for k, x in enumerate(side) if
                                    pymbt.analysis.tm(seqs[k]) > tm]
                else:
                    seqs = [targets[i][-x:] for x in side]
                    report[i][j] = [x for k, x in enumerate(side) if
                                    pymbt.analysis.tm(seqs[k]) > tm]

        # 3. See if there's more than one result. If so, throw exception
        count = len([z for x in report for y in x for z in y])
        if not count:
            raise ValueError('Failed to find compatible Gibson ends.')
        elif count > 1:
            raise AmbiguousGibsonError('multiple compatible ends.')

        # 4. There must be one result. Where is it?
        for i, target in enumerate(report):
            for j, side in enumerate(target):
                if side:
                    match = (i, j, side[0])

        # 5. Combine pieces together
        left_side = self._working_list[0]
        right_side = self._working_list.pop(match[0] + 1)
        if match[1] == 1:
            right_side = right_side.reverse_complement()

        self._working_list[0] = left_side + right_side[match[2]:]

    def _fuse_last(self, homology, tm):
        '''With one sequence left, attempt to fuse it to itself.

        :param homology: length of terminal homology in bp.
        :type homology: int
        :raises: AmbiguousGibsonError if either of the termini are palindromic
                 (would bind self-self).
                 ValueError if the ends are not compatible.

        '''
        # 1. Get report on self-self
        pattern = self._working_list[0]
        report = homology_report(pattern, pattern)
        report_l = [x for x in report[0] if x != len(pattern)]
        report_l = [x for x in report_l if x > homology]
        report_l = [x for x in report_l if
                    pymbt.analysis.tm(pattern[:x + 1]) > tm]
        if not report_l:
            raise ValueError('Failed to find compatible Gibson ends.')
        elif len(report_l) > 1:
            raise AmbiguousGibsonError('multiple compatible ends.')
        self._working_list[0] = pattern[:-report_l[0]]
        self._working_list[0] = self._working_list[0].circularize()


def homology_report(seq1, seq2):
    '''Given two sequence inputs, find 3\' revcomp identities.

    :param seq1: Sequence to test 3\' end of Watson strand of
    :type seq1: pymbt.DNA
    :param seq2: Sequence to test both strands of, 3\' end
    :type seq1: pymbt.DNA
    :returns: List of left and right identities.
    :rtype: list of ints

    '''
    # Go through each other sequence, both forward and revcomp sequences,
    # and check whether first base matches last / first of initial sequence
    reference = str(seq1)
    watson_3 = str(seq2)
    crick_3 = str(seq2.reverse_complement())

    left_list = []
    for i in range(min(len(reference), len(watson_3))):
        if reference[-(i + 1):] == watson_3[:i + 1]:
            left_list.append((i + 1))

    right_list = []
    for i in range(min(len(reference), len(crick_3))):
        if reference[-(i + 1):] == crick_3[:i + 1]:
            right_list.append((i + 1))

    return [left_list, right_list]
