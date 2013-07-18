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
    # Find match in top or bottom strands for each primere
    p1_matches = template.locate(primer1.anneal)
    p2_matches = template.locate(primer2.anneal)
    # Combine the matches
    combined = [p1 + p2 for p1, p2 in zip(p1_matches, p2_matches)]
    # Make sure there's no ambiguities
    for strand in combined:
        if len(strand) > 1:
            raise Exception('Ambiguous priming detected.')
    # Make 'reverse' index useful for slicing
    combined[1][0] = len(template) - combined[1][0]
    # Subset
    return template[combined[0][0]:combined[1][0]]
