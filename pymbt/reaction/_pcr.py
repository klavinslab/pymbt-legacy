'''PCR reaction(s).'''


def pcr(template, primer1, primer2):
    '''Simulate a PCR (no support for ambiguous PCRs).

    :param template: DNA template from which to PCR.
    :type template: pymbt.sequence.DNA
    :param primer1: First PCR primer.
    :type primer1: pymbt.sequence.Primer
    :param primer2: First PCR primer.
    :type primer2: pymbt.sequence.Primer

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
    # if rev == fwd (start/stop is the same), amplify whole plasmid
    # if rev > fwd and rev - fwd is less than the size of one of the primers,
    # then a weird overlap starts to occur - primer dimers become likely.
    # if rev > fwd, rev - fwd < maxlen(fwd, rev), weird stuff happens:
    #   if rev - fwd >= minlen(fwd, rev) as well, primer dimers are guaranteed
    #   if rev - fwd < minlen(fwd, rev), primer dimers are also likely, but
    # less so.
    # The problem of primer dimers is separate from just trying to simulate a
    # PCR based on molecular genetics rules. What's the best compromise to use
    # here?
    #   If rev - fwd > maxlen(rev, fwd), everything is OK.
    #   If rev - fwd == maxlen(rev, fwd), should ampify primer dimer (desired?)
    #   If rev - fwd < maxlen(rev, fwd), there's a risk of weird amplifications
    #       Should raise exception for now with a note to address it later
    # Subset
    if rev - fwd < max(len(primer1), len(primer2)):
        raise Exception('Primers bind in one another, no solution implemented')
    if rev > fwd:
        amplicon = template[fwd:rev]
    else:
        amplicon = template[fwd:] + template[:rev]
    if primer1.overhang:
        amplicon = primer1.overhang.set_stranded('ds') + amplicon
    if primer2.overhang:
        amplicon += primer2.overhang.set_stranded('ds').reverse_complement()
    return amplicon
