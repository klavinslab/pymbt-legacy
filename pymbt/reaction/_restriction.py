'''
Restriction endonuclease reactions.

'''


class Restriction(object):
    '''
    Restriction endonuclease reaction.

    '''

    def __init__(self, dna, restriction_site):
        '''
        :param dna: DNA template to digest.
        :type dna: pymbt.sequence.DNA
        :param restriction_site: Restriction site to use.
        :type restriction_site: RestrictionSite

        '''

        self.template = dna
        self.restriction_site = restriction_site
        self.output = []

    def run(self):
        '''
        Execute Restriction and return results.

        '''

        # TODO: search for top *and* bottom strand match - currently succeeds
        # if only one matches, but enzymes don't work that way
        # TODO: Figure out what to do with ambiguities.
        #       example: site with both AsuI (TTCGAA) and EcoRI (GAATTC) sites
        #       if EcoRI cuts first the AsuI site is gone, and vice versa
        #       basically need to enumerate all different orders in which
        #       enzymes could cut... which is potentially a lot...
        #           product ( site_n! ) where n is the number of times
        #           that the site appears in the sequence
        #           So if EcoRI cuts 3 times and AsuI cuts twice, there are
        #           3! * 2! = 12 permutations
        #           if EcoRI cuts 5 times and AsuI cuts 4 times, there are
        #           5! * 4! = 120 * 24 = 2880 permutations.
        #           Hopefully these operations are very fast!
        dna_pattern = str(self.restriction_site.recognition_site)
        both_indices = self.template.locate(dna_pattern)
        if both_indices:
            if both_indices[0]:
                cut1 = both_indices[0][0]
                return self._cut(cut1)
        else:
            print 'Restriction site not present in template DNA.'
            return [self.template]

        # Check for circular - if so cut once then proceed

        # Split into new pieces including full restriction site on cut end

        # Resect correct amount to generate sticky ends

    def _cut(self, index):
        '''
        Cuts template once at the specified index.

        '''

        dna_list = []

        cut_site = self.restriction_site.cut_site
        leftmost = index + min(cut_site)

        top_cut = leftmost + cut_site[0]
        bottom_cut = leftmost + cut_site[1]

        if self.template.type == 'circular':
            # List contains only linearized DNA
            dna_list.append(self.template.linearize(leftmost))
        else:
            # List has left piece, then right piece in digest
            dna_list.append(self.template[0:leftmost])
            dna_list.append(self.template[leftmost:])

        if top_cut != bottom_cut:
            # Find sequence between cut sites
            between = self.template[leftmost:index + max(cut_site)]
            between_rev = between.reverse_complement()

            # Prepare ssDNA versions of overhangs
            top_overhang = between.five_resect()
            bottom_overhang_rev = between_rev.five_resect()
            bottom_overhang = bottom_overhang_rev.reverse_complement()

            # If cut has 3' overhang, reverse complement overhangs
            if top_cut > bottom_cut:
                top_overhang = top_overhang.reverse_complement()
                bottom_overhang = bottom_overhang.reverse_complement()
            if self.template.type == 'circular':
                dna_list[0] = dna_list[0][len(between):]
                dna_list[0] += top_overhang
                dna_list[0] = bottom_overhang + dna_list[0]
            else:
                dna_list[0] += top_overhang
                dna_list[1] = dna_list[1][len(between):]
                dna_list[1] = bottom_overhang + dna_list[1]

        return dna_list


# Inverted repeat checking - should be done here, not core DNA class
#        inverted_repeat = analysis.palindrome(pattern)
#        if inverted_repeat:
#            # subtract all occurrences in top from bottom
#            subtract = [len(self.top) - index - len(pattern) for index in
#                        top_starts]
#            bottom_starts = [x for x in subtract if x not in bottom_starts]
