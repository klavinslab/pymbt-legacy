'''PCR reaction(s).'''


class PCR(object):
    '''Simulate a PCR (no support for ambiguous PCRs).'''
    def __init__(self, primer1, primer2, template):
        self.primer1 = primer1
        self.primer2 = primer2
        self.template = template

    def run(self):
        # Find match in top or bottom strands for each primere
        p1_matches = self.template.locate(self.primer1.anneal)
        p2_matches = self.template.locate(self.primer2.anneal)
        # Combine the matches
        combined = [p1 + p2 for p1, p2 in zip(p1_matches, p2_matches)]
        # Make sure there's no ambiguities
        for strand in combined:
            if len(strand) > 1:
                raise Exception('Ambiguous priming detected.')
        # Make 'reverse' index useful for slicing
        combined[1][0] = len(self.template) - combined[1][0]
        # Subset
        return self.template[combined[0][0]:combined[1][0]]
