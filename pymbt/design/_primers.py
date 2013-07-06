'''
Primer design tools.

'''

from pymbt import analysis


# TODO: combine DesignPrimer and DesignPrimer gene in a sane way.
# ideas:
#   Give DesignPrimer a 'reverse' option that computes a reverse primer
#   Give DesignPrimer a 'both' option that computes both
class DesignPrimer(object):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

    '''

    def __init__(self, dna_object, tm=72, min_len=10, tm_undershoot=1,
                 tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                 overhang=None):
        '''
        :param dna_object: object of type DNA.
        :type dna_object: DNA
        :param tm: Ideal primer Tm in degrees C.
        :type tm: float
        :param min_len: Minimum primer length.
        :type min_len: int
        :param tm_undershoot: Allowed Tm undershoot.
        :type tm_undershoot: float
        :param tm_overshoot: Allowed Tm overshoot.
        :type tm_overshoot: float
        :param end_gc: Obey the 'end on G or C' rule.
        :type end_gc: bool
        :param tm_method: Melting temp calculator method to use.
        :type tm_method: string
        :param overhang: Append the primer to this overhang sequence.
        :type overhang: DNA

        '''

        # TODO: deal with sticky ended inputs or require a new DNA type
        # that can't have them (dsDNA)

        self.template = dna_object

        # Storing inputs as attributes so they can be modified in-place?
        self.tm = tm
        self.min_len = min_len
        self.tm_undershoot = tm_undershoot
        self.tm_overshoot = tm_overshoot
        self.end_gc = end_gc
        self.tm_method = tm_method
        self.overhang = overhang

    def run(self):
        '''
        Execute the design algorithm.

        '''
        primer, tm = _design_primer(self.template, self.tm, self.min_len,
                                    self.tm_undershoot, self.tm_overshoot,
                                    self.end_gc, self.tm_method, self.overhang)
        return primer, tm


class DesignPrimerGene(object):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

    '''

    def __init__(self, dna_object, tm=72, min_len=10, tm_undershoot=1,
                 tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                 overhangs=None):
        '''
        :param dna_object: object of type DNA.
        :type dna_object: DNA
        :param tm: Ideal primer Tm in degrees C.
        :type tm: float
        :param min_len: Minimum primer length.
        :type min_len: int
        :param tm_undershoot: Allowed Tm undershoot.
        :type tm_undershoot: float
        :param tm_overshoot: Allowed Tm overshoot.
        :type tm_overshoot: float
        :param end_gc: Obey the 'end on G or C' rule.
        :type end_gc: bool
        :param tm_method: Melting temp calculator method to use.
        :type tm_method: string
        :param overhangs: Sequences to append to the primers (2-tuple).
        :type overhangs: tuple of DNA object

        '''

        # TODO: deal with sticky ended inputs or require a new DNA type
        # that can't have them (dsDNA)

        # Type checking on input
        self.template = dna_object

        # Storing inputs as attributes so they can be modified in-place?
        self.tm = tm
        self.min_len = min_len
        self.tm_undershoot = tm_undershoot
        self.tm_overshoot = tm_overshoot
        self.end_gc = end_gc
        self.tm_method = tm_method
        self.overhangs = overhangs

    def run(self):
        '''
        Execute the design algorithm.

        '''
        template = self.template
        primers_list = _design_primer_gene(template, self.tm, self.min_len,
                                           self.tm_undershoot,
                                           self.tm_overshoot, self.end_gc,
                                           self.tm_method, self.overhangs)
        return primers_list


def _design_primer(dna_object, tm=72, min_len=10, tm_undershoot=1,
                   tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                   overhang=None):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

    :param dna_object: Input DNA.
    :type dna_object: DNA object.
    :param tm: Ideal primer Tm in degrees C.
    :type tm: float
    :param min_len: Minimum primer length.
    :type min_len: int
    :param tm_undershoot: Allowed Tm undershoot.
    :type tm_undershoot: float
    :param tm_overshoot: Allowed Tm overshoot.
    :type tm_overshoot: float
    :param end_gc: Obey the 'end on G or C' rule.
    :type end_gc: bool
    :param tm_method: Melting temp calculator method to use.
    :type tm_method: string
    :param overhang: Append the primer to this overhang sequence.
    :type overhang: str

    '''

    # Check Tm of input sequence to see if it's already too low
    seq_tm = analysis.tm(dna_object, method=tm_method)
    if seq_tm < (tm - tm_undershoot):
        msg = 'Input sequence Tm is lower than primer Tm setting'
        raise Exception(msg)

    # Focus on first 90 bases - shouldn't need more than 90bp to anneal
    dna_object = dna_object[0:90]

    # Generate primers from min_len to 'tm' + tm_overshoot
    primers_tms = []

    last_tm = 0
    max_tm = tm + tm_overshoot
    bases = min_len

    while last_tm <= max_tm and (bases != len(dna_object)):
        new_primer = dna_object[0:bases]
        last_tm = analysis.tm(new_primer, method=tm_method)
        primer_tm = (new_primer, last_tm)
        primers_tms.append(primer_tm)
        bases += 1

    # Trim primer list based on tm_undershoot and end_gc
    tmin = tm - tm_undershoot

    primers = []
    tms = []
    primers_tms = [(primer, melt) for primer, melt in primers_tms if
                   melt >= tmin]

    if end_gc:
        primers_tms = [(primer, melt) for primer, melt in primers_tms if
                       primer.endswith(('c', 'g'))]

    primers, tms = zip(*primers_tms)

    if not primers:
        raise Exception('No primers could be generated using these settings')

    # Find the primer closest to the set Tm
    tm_diffs = [abs(x - tm) for x in tms]
    best_index = tm_diffs.index(min(tm_diffs))
    best_primer = primers[best_index]
    best_tm = tms[best_index]

    # Make it single-stranded
    best_primer = best_primer.set_stranded('ss')

    if overhang:
        overhang = overhang.set_stranded('ss')
        best_primer = overhang + best_primer

    return best_primer, best_tm


def _design_primer_gene(dna_object, tm=72, min_len=10, tm_undershoot=1,
                        tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                        overhangs=None):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

    :param seq: Input sequence.
    :type seq: str
    :param tm: Ideal primer Tm in degrees C.
    :type tm: float
    :param min_len: Minimum primer length.
    :type min_len: int
    :param tm_undershoot: Allowed Tm undershoot.
    :type tm_undershoot: float
    :param tm_overshoot: Allowed Tm overshoot.
    :type tm_overshoot: float
    :param end_gc: Obey the 'end on G or C' rule.
    :type end_gc: bool
    :param tm_method: Melting temp calculator method to use.
    :type tm_method: string
    :param overhangs: 2-tuple of overhang sequences.
    :type overhangs: tuple

    '''
    if not overhangs:
        overhangs = [None, None]

    templates = [dna_object, dna_object.reverse_complement()]
    primer_list = []

    for template, overhang in zip(templates, overhangs):
        primer, tm = _design_primer(template, tm=tm, min_len=min_len,
                                    tm_undershoot=tm_undershoot,
                                    tm_overshoot=tm_overshoot, end_gc=end_gc,
                                    tm_method=tm_method,
                                    overhang=overhang)
        primer_list.append((primer, tm))
    return primer_list