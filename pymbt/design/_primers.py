'''Primer design tools.'''

from pymbt import analysis
from pymbt import sequence
import warnings


# TODO: combine DesignPrimer and DesignPrimerGene
# ideas:
#   Give DesignPrimer a 'reverse' option that computes a reverse primer
#   Give DesignPrimer a 'both' option that computes both
#   Provide DesignPrimer with 'forward' and 'reverse' methods
class DesignPrimer(object):
    '''Design primer to a nearest-neighbor Tm setpoint.'''
    def __init__(self, dna, tm=72, min_len=10, tm_undershoot=1,
                 tm_overshoot=3, end_gc=False, tm_parameters='cloning',
                 overhang=None):
        '''
        :param dna: Sequence for which primers will be designed.
        :type dna: pymbt.sequence.DNA
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
        :param tm_parameters: Melting temp calculator method to use.
        :type tm_parameters: string
        :param overhang: Append the primer to this overhang sequence.
        :type overhang: DNA

        '''
        # TODO: deal with sticky ended inputs or require a new DNA type
        # that can't have them (dsDNA)

        self.template = dna

        # Storing inputs as attributes so they can be modified in-place?
        self.tm = tm
        self.min_len = min_len
        self.tm_undershoot = tm_undershoot
        self.tm_overshoot = tm_overshoot
        self.end_gc = end_gc
        self.tm_parameters = tm_parameters
        self.overhang = overhang

        self.anneal = None
        self.primer = None
        self.primer_tm = None
        self.high_structure = False

    def run(self, check_structure=True):
        '''Design the primer.'''
        self.primer = _design_primer(self.template, self.tm,
                                     self.min_len,
                                     self.tm_undershoot,
                                     self.tm_overshoot,
                                     self.end_gc,
                                     self.tm_parameters,
                                     self.overhang)
        if check_structure:
            self.structure()
        return self.primer

    def structure(self):
        '''Check annealing sequence for structure.'''
        # Check whole primer for high-probability structure, focus in on
        # annealing sequence, report average
        nupack = analysis.Nupack(self.primer.primer)
        pairs = nupack.pairs()
        anneal_len = len(self.primer.anneal)
        pairs_mean = sum(pairs[-anneal_len:]) / anneal_len
        if pairs_mean < 0.5:
            self.high_structure = True
            warnings.warn('High probability structure', Warning)
        return pairs_mean


class DesignPrimerGene(object):
    '''Design two primers for cloning a gene (or any arbitrary sequence).'''
    def __init__(self, dna, tm=72, min_len=10, tm_undershoot=1,
                 tm_overshoot=3, end_gc=False, tm_parameters='cloning',
                 overhangs=None):
        '''
        :param dna: Sequence for which to design primers.
        :type dna: pymbt.sequence.DNA
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
        :param tm_parameters: Melting temp calculator method to use.
        :type tm_parameters: string
        :param overhangs: Sequences to append to the primers (2-tuple).
        :type overhangs: tuple of DNA object

        '''
        # TODO: deal with sticky ended inputs or require a new DNA type
        # that can't have them (dsDNA)

        # Type checking on input
        self.template = dna

        # Storing inputs as attributes so they can be modified in-place?
        self.tm = tm
        self.min_len = min_len
        self.tm_undershoot = tm_undershoot
        self.tm_overshoot = tm_overshoot
        self.end_gc = end_gc
        self.tm_parameters = tm_parameters
        self.overhangs = overhangs

    def run(self):
        '''Design the primers.'''
        template = self.template
        primers_list = _design_primer_gene(template, self.tm, self.min_len,
                                           self.tm_undershoot,
                                           self.tm_overshoot, self.end_gc,
                                           self.tm_parameters, self.overhangs)
        return primers_list


def _design_primer(dna, tm=72, min_len=10, tm_undershoot=1,
                   tm_overshoot=3, end_gc=False, tm_parameters='cloning',
                   overhang=None):
    '''Design primer to a nearest-neighbor Tm setpoint.

    :param dna: Sequence for which to design a primer.
    :type dna: pymbt.sequence.DNA
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
    :param tm_parameters: Melting temp calculator method to use.
    :type tm_parameters: string
    :param overhang: Append the primer to this overhang sequence.
    :type overhang: str

    '''
    # Check Tm of input sequence to see if it's already too low
    seq_tm = analysis.tm(dna, parameters=tm_parameters)
    if seq_tm < (tm - tm_undershoot):
        msg = 'Input sequence Tm is lower than primer Tm setting'
        raise ValueError(msg)
    # Focus on first 90 bases - shouldn't need more than 90bp to anneal
    dna = dna[0:90]

    # Generate primers from min_len to 'tm' + tm_overshoot
    # TODO: this is a good place for optimization. Only calculate as many
    # primers as are needed. Use binary search.
    primers_tms = []
    last_tm = 0
    bases = min_len
    while last_tm <= tm + tm_overshoot and bases != len(dna):
        next_primer = dna[0:bases]
        last_tm = analysis.tm(next_primer, parameters=tm_parameters)
        primers_tms.append((next_primer, last_tm))
        bases += 1

    # Trim primer list based on tm_undershoot and end_gc
    primers_tms = [(primer, melt) for primer, melt in primers_tms if
                   melt >= tm - tm_undershoot]
    if end_gc:
        primers_tms = [(primer, melt) for primer, melt in primers_tms if
                       primer.endswith(('c', 'g'))]
    if not primers_tms[0]:
        raise Exception('No primers could be generated using these settings')

    # Find the primer closest to the set Tm, make it single stranded
    tm_diffs = [abs(melt - tm) for primer, melt in primers_tms]
    best_index = tm_diffs.index(min(tm_diffs))
    best_primer, best_tm = primers_tms[best_index]
    best_primer = best_primer.set_stranded('ss')

    # Apply overhang
    if overhang:
        overhang = overhang.set_stranded('ss')
        best_primer = overhang + best_primer

    return sequence.Primer(best_primer, overhang, best_tm)


def _design_primer_gene(dna, tm=72, min_len=10, tm_undershoot=1,
                        tm_overshoot=3, end_gc=False,
                        tm_parameters='cloning',
                        overhangs=None):
    '''Design primers for cloning any arbitrary sequence..

    :param dna: Input sequence.
    :type dna: pymbt.sequence.DNA
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
    :param tm_parameters: Melting temp calculator method to use.
    :type tm_parameters: string
    :param overhangs: 2-tuple of overhang sequences.
    :type overhangs: tuple

    '''
    if not overhangs:
        overhangs = [None, None]
    templates = [dna, dna.reverse_complement()]
    primer_list = []
    for template, overhang in zip(templates, overhangs):
        primer = _design_primer(template, tm=tm, min_len=min_len,
                                tm_undershoot=tm_undershoot,
                                tm_overshoot=tm_overshoot, end_gc=end_gc,
                                tm_parameters=tm_parameters,
                                overhang=overhang)
        primer_list.append(primer)
    return primer_list
