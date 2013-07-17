'''Gibson reaction simulation.'''
# TODO: allow production of a linear fragment as well


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

    def run(self, homology_min=15):
        '''Run the Gibson reaction.'''
        # TODO: could miss an ambiguity that has even more homology to a given
        # piece than the first one that matches perfectly (index==0)

        # Until there's only one piece left, start fusing them together
        while len(self._seq_list) > 1:
            self._find_fuse_next(homology_min)
        # Fuse the last piece to itself
        self._fuse_last(homology_min)
        return self._seq_list[0]

    def _find_fuse_next(self, homology):
        print 'adding'
        while True:
            current = self._seq_list[0]
            to_locate = current.top[-homology:]
            found = []
            # Generate matches for current to_locate sequence
            for i, x in enumerate(self._seq_list):
                if x != current:
                    found.append(x.locate(to_locate))
            print found
            # If there are no results, throw an exception
            flat_results = [z for x in found for y in x for z in y]
            if not flat_results:
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
            # If more than one are perfect matches, Gibson is ambiguous
            if matches:
                if len(matches) > 1:
                    msg = 'Ambiguous Gibson (multiple compatible ends)'
                    raise Exception(msg)
                else:
                    # If all those other checks passed, match is unique!
                    # Combine the current sequence with that sequence
                    matched = self._seq_list.pop(matches[0][0] + 1)
                    # Was it top or bottom?
                    if matches[0][1] == 1:
                        print 'revcomping'
                        matched = matched.reverse_complement()
                    print 'current template len: {}'.format(len(current))
                    print 'next piece len: {}'.format(len(matched))
                    new_current = current[0:-homology] + matched
                    self._seq_list[0] = new_current
                    break
            homology += 1

    def _fuse_last(self, homology):
        # Should use basically the same approach as above...
        print 'completing'
        while True:
            current = self._seq_list[0]
            to_locate = current.top[-homology:]
            # Generate matches for current to_locate sequence
            found = current.locate(to_locate)[0]
            print found
            if not found:
                raise Exception('Failed to find compatible Gibson ends.')
            # There are matches so see if any match homology length
            # If more than one are perfect matches, Gibson is ambiguous
            matches = []
            for index in found:
                if index == 0:
                    matches.append(index)
            if len(matches) > 1:
                msg = 'Ambiguous Gibson (multiple compatible ends)'
                raise Exception(msg)
            elif len(matches) == 1:
                # If all those other checks passed, match is unique!
                # Trim off homology length and circularize
                new_current = current[0:-homology]
                self._seq_list[0] = new_current
                self._seq_list[0] = self._seq_list[0].circularize()
                break
            homology += 1
