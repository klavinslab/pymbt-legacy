"""Simulate building a construct by assembling oligos with PCA."""


class AssemblyError(Exception):
    """Raise exception if assembly can't complete for any reason."""
    pass


def assemble_oligos(dna_list):
    """Given a list of DNA sequences, assemble into a single construct.
    :param dna_list: List of DNA sequences - they must be single-stranded.
    :type dna_list: pymbt.sequence.DNA list
    :raises: AssemblyError if it can't assemble for any reason.
    :returns: A single assembled DNA sequence
    :rtype: pymbt.sequence.DNA

    """
    # Find all matches for every oligo. If more than 2 per side, error.
    # Self-oligo is included in case the 3' end is self-complementary.
    # 3' end matches (reference is 3' end, query is 5' end of other oligo)
    match_3 = [match(seq, dna_list, right=True) for i, seq in
               enumerate(dna_list)]
    # 5' end matches (reference is 5' end, query is 3' end of other oligo)
    match_5 = [match(seq, dna_list, right=False) for i, seq in
               enumerate(dna_list)]
    # Assemble into 2-tuple
    zipped = zip(match_3, match_5)
    return zipped
    for i, oligo_matches in enumerate(zipped):
        if not any(oligo_matches):
            error = "Oligo {} doesn't bind any others.".format(i + 1)
            raise AssemblyError(error)
    # Isolate oligos with only one match - there should be exactly 2
    # < 2 implies missing oligo. >2 implies an extra oligo
    # Treat first one found (left to right from input list) as 'first' oligo.
    ends = [i for i, x in enumerate(zipped) if any([y is None for y in x])]
    if len(ends) > 2:
        raise AssemblyError("More than 2 end oligos found.")
    if len(ends) < 2:
        raise AssemblyError("Fewer than 2 end oligos found.")
    return ends
    # Chain first oligo to second oligo, fill in side that won't be used
    # as overhang.
    # If a match can't be found
    # Repeat until assembly is complete.


def match(reference, query, right=True, start_size=12):
    """Given reference oligo, find queries that uniquely match its 3' end"""
    size = start_size
    found = []
    # Reverse complementing here provides massive speedup
    rev_query = [seq.reverse_complement().flip() for seq in query]
    while not found and not size > len(reference):
        for i, seq in enumerate(rev_query):
            if right:
                if reference[-size:] == seq[:size]:
                    found.append(i)
            else:
                if reference[:size] == seq[-size:]:
                    found.append(i)
        size += 1
    if len(found) > 1:
        raise AssemblyError("Ambiguous oligo binding")
    if not found:
        return None
    else:
        return found
