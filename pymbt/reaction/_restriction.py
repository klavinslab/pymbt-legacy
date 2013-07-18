'''Restriction endonuclease reactions.'''


# FIXME: pmod4G-yevenus-stuff is incorrect when digested with NcoI
def digest(dna, restriction_enzyme):
    '''Restriction endonuclease reaction.

    :param dna: DNA template to digest.
    :type dna: pymbt.sequence.DNA
    :param restriction_site: Restriction site to use.
    :type restriction_site: RestrictionSite

    '''
    pattern = restriction_enzyme.recognition_site
    located = dna.locate(pattern)
    if not located[0] and not located[1]:
        # TODO: should this raise an exception?
        print 'Restriction site not present in template DNA.'
        return [dna]
    # Bottom strand indices are relative to the bottom strand 5' end.
    # Convert to same type as top strand
    pattern_len = len(pattern)
    r_indices = [len(dna) - index - pattern_len for index in
                 located[1]]
    # If sequence is palindrome, remove redundant results
    if pattern.is_palindrome():
        r_indices = [index for index in r_indices if index not in
                     located[0]]
    # Flatten cut site indices
    cut_sites = sorted(located[0] + r_indices)
    # Go through each cut site starting at highest one
    # Cut remaining template once, generating remaining + new
    current = [dna]
    for cut_site in cut_sites[::-1]:
        new = _cut(current, cut_site, restriction_enzyme)
        current.append(new[1])
        current.append(new[0])
    current.reverse()
    # Combine first and last back together if digest was circular
    if dna.topology == 'circular':
        current[0] = current.pop() + current[0]
    return current


def _cut(dna, index, restriction_enzyme):
    '''Cuts template once at the specified index.

    :param dna: DNA to cut
    :type dna: pymbt.sequence.DNA
    :param index: index at which to cut
    :type index: int
    :param restriction_enzyme: Enzyme with which to cut
    :type restriction_enzyme: pymbt.sequence.RestrictionSite

    '''
    # Find absolute indices at which to cut
    cut_site = restriction_enzyme.cut_site
    top_cut = index + cut_site[0]
    bottom_cut = index + cut_site[1]

    # Isolate left and ride sequences
    to_cut = dna.pop()
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
        left = left.three_resect(abs(diff))
        right_r = right.reverse_complement()
        right = right_r.three_resect(abs(diff)).reverse_complement()

    return [left, right]
