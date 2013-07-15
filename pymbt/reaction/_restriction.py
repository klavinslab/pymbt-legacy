'''Restriction endonuclease reactions.'''

from pymbt.analysis import palindrome


class Restriction(object):
    '''Restriction endonuclease reaction.'''

    def __init__(self, dna, restriction_site):
        '''
        :param dna: DNA template to digest.
        :type dna: pymbt.sequence.DNA
        :param restriction_site: Restriction site to use.
        :type restriction_site: RestrictionSite

        '''
        self.template = dna
        self.current = []
        self.restriction_site = restriction_site
        self.output = []
        self.ran = False

    def run(self):
        '''Simulate the restriction digest and return results.'''
        pattern = self.restriction_site.recognition_site
        located = self.template.locate(pattern)
        if not located[0] and not located[1]:
            print 'Restriction site not present in template DNA.'
            return [self.template]
        # Bottom strand indices are relative to the bottom strand 5' end.
        # Convert to same type as top strand
        pattern_len = len(pattern)
        r_indices = [len(self.template) - index - pattern_len for index in
                     located[1]]
        # If sequence is palindrome, remove redundant results
        if palindrome(pattern):
            r_indices = [index for index in r_indices if index not in
                         located[0]]
        # Flatten cut site indices
        cut_sites = sorted(located[0] + r_indices)
        # Go through each cut site starting at highest one
        # Cut remaining template once, generating remaining + new
        self.current = [self.template]
        for cut_site in cut_sites[::-1]:
            new = self._cut(cut_site)
            self.current.append(new[1])
            self.current.append(new[0])
        self.current.reverse()
        # Combine first and last back together if digest was circular
        if self.template.topology == 'circular':
            print 'yep'
            self.current[0] += self.current.pop()
        return self.current

    def _cut(self, index):
        '''Cuts template once at the specified index.

        :param index: index at which to cut
        :type index: int

        '''
        # Find absolute indices at which to cut
        cut_site = self.restriction_site.cut_site
        top_cut = index + cut_site[0]
        bottom_cut = index + cut_site[1]

        # Isolate left and ride sequences
        to_cut = self.current.pop()
        left = to_cut[0:top_cut]
        right = to_cut[bottom_cut:]

        # If applicable, leave overhangs
        diff = top_cut - bottom_cut
        if not diff:
            # Blunt-end cutter, no adjustment necessaryy
            pass
        elif diff > 0:
            # 3' overhangs
            left_r = left.reverse_complement()
            left = left_r.five_resect(diff).reverse_complement()
            right = right.five_resect(diff)
        else:
            # 5' overhangs
            left = left.three_resect(diff)
            right_r = right.reverse_complement()
            right = right_r.three_resect(diff).reverse_complement()

        return [left, right]
