'''PCR reaction(s).'''


def pcr(template, primer1, primer2):
    '''Simulate a PCR (no support for ambiguous PCRs).

    :param template: DNA template from which to PCR.
    :type template: pymbt.sequence.DNA
    :param primer1: First PCR primer.
    :type primer1: pymbt.sequence.Primer
    :param primer2: First PCR primer.
    :type primer2: pymbt.sequence.Primer
    :returns: A dsDNA Amplicon.
    :rtype: pymbt.sequence.DNA
    :raises: Exception if a primer binds more than once on the template.
             Exception if primers bind in overlapping sequence of the template.
             Exception if the PCR would work on a circular version of the
             template (implies that input was linear).

    '''
    # Find match in top or bottom strands for each primer
    p1_matches = template.locate(primer1.anneal)
    p2_matches = template.locate(primer2.anneal)
    fwds = p1_matches[0] + p2_matches[0]
    revs = p2_matches[1] + p1_matches[1]
    # Make sure there's no ambiguities
    if len(fwds) > 1 or len(revs) > 1:
        raise Exception('Ambiguous priming detected.')
    # Make 'reverse' index useful for slicing
    fwd = fwds[0]
    rev = len(template) - revs[0]
    # TODO: circular search will muck things up. If primers are at the very
    # beginning or end of the plasmid coordinates, things will get weird
    # TODO: Should actually just evaluate the length of the product prior
    # to adding primer overhangs, compare to length of anneal regions.
    #   But would need to keep track of sign - fwd - rev != rev - fwd
    '''Notes on PCR amplification decisions:
        If rev == fwd, primers should ampify entire plasmid
        If rev - fwd >= max(len(rev), len(fwd)), amplify sequence normally
        If rev - fwd < max(len(rev), len(fwd)), raise exception - who knows
        how this construct will amplify

    '''
    if rev - fwd < max(len(primer1), len(primer2)) and rev - fwd > 0:
        raise Exception('Primers bind in one another, no solution implemented')

    # Subset
    if rev > fwd:
        amplicon = template[fwd:rev]
    else:
        if template.topology == 'circular':
            amplicon = template[fwd:] + template[:rev]
        else:
            raise Exception('Primers would amplify if template were circular.')
    if primer1.overhang:
        amplicon = primer1.overhang.set_stranded('ds') + amplicon
    if primer2.overhang:
        amplicon += primer2.overhang.set_stranded('ds').reverse_complement()
    return amplicon
