'''Gibson reaction simulation.'''


class Gibson(object):
    '''Gibson reaction.'''
    def __init__(self, seq_list):
        '''
        :param seq_list: list of DNA sequences to Gibson
        :type seq_list: list of pymbt.sequence.DNA

        '''
        self._seq_list = seq_list
        if any(seq.topology == 'circular' for seq in self._seq_list):
            raise ValueError('Input sequences must be linear, not circular.')

    def run(self, linear=False, homology_min=15):
        '''Run the Gibson reaction.

        :param linear: Return a linear fragment rather than attempting to
                       generate a circular construct.
        :type linear: bool
        :param homology_min: minimum bp of homology allowed
        :type homology_min: int

        '''
        # FIXME: could miss an ambiguity that has even more homology to a given
        # piece than the first one that matches perfectly (index==0)
        # Strategy: Once matches are found on other strands, start expanding
        # them until the end of the sequence is reached. If *that* matches
        # the original sequence homology area (but adjusted to be larger),
        # keep it if there's only one, if there's more raise error

        # Attempt to fuse fragments together until only one is left
        while len(self._seq_list) > 1:
            self._find_fuse_next(homology_min)
        if not linear:
            # Fuse the final fragment to itself
            self._fuse_last(homology_min)
        return self._seq_list[0]

    def _find_fuse_next(self, homology):
        while True:
            pattern = self._seq_list[0].top[-homology:]
            # Generate matches for working sequence terminal homology
            found = [x.locate(pattern) for x in self._seq_list if
                     x != self._seq_list[0]]
            # If there are no results, throw an exception
            if not [z for x in found for y in x for z in y]:
                raise Exception('Failed to find compatible Gibson ends.')
            # If there are results, see if any match homology length
            matches = []
            for i, (top, bottom) in enumerate(found):
                for index in top:
                    if index == 0:
                        matches.append((i, 0))
                for b_index in bottom:
                    if b_index == 0:
                        matches.append((i, 1))
            if matches:
                # If more than one are perfect matches, Gibson is ambiguous
                if len(matches) > 1:
                    msg = 'Ambiguous Gibson (multiple compatible ends)'
                    raise Exception(msg)
                else:
                    # If all those other checks passed, match is unique!
                    # Combine the working sequence with that sequence
                    matched = self._seq_list.pop(matches[0][0] + 1)
                    # Reverse complement before concatenation if bottom strand
                    if matches[0][1] == 1:
                        matched = matched.reverse_complement()
                    self._seq_list[0] = self._seq_list[0][:-homology] + matched
                    break
            homology += 1

    def _fuse_last(self, homology):
        # Should use basically the same approach as above...
        while True:
            pattern = self._seq_list[0].top[-homology:]
            # Generate matches for working sequence terminal homology
            found = self._seq_list[0].locate(pattern)[0]
            found.pop(-1)
            if not found:
                raise Exception('Failed to find compatible Gibson ends.')
            # There are matches so see if any match homology length
            # If more than one are perfect matches, Gibson is ambiguous
            matches = [index for index in found if index == 0]
            if len(matches) > 1:
                msg = 'Ambiguous Gibson (multiple compatible ends)'
                raise Exception(msg)
            elif len(matches) == 1:
                # If all those other checks passed, match is unique!
                # Trim off redundant homology and circularize
                self._seq_list[0] = self._seq_list[0][:-homology].circularize()
                break
            homology += 1