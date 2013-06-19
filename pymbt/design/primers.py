'''Primer design tools.'''

from pymbt import sequence
from pymbt import analysis
from pymbt.sequence.utils import check_alphabet
from pymbt.sequence.utils import reverse_complement


class DesignPrimer(object):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

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
    :type overhang: str

    '''

    def __init__(self, dna_object, tm=72, min_len=10, tm_undershoot=1,
                 tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                 overhang=''):
        # TODO: type checking
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
        primer = _design_primer(str(self.template), self.tm, self.min_len,
                                self.tm_undershoot, self.tm_overshoot,
                                self.end_gc, self.tm_method, self.overhang)
        #primer = _design_primer(str(self.template), **self.kwargs)
        return sequence.DNA(primer[0], stranded='ss')


class DesignPrimerGene(object):
    '''
    Design primer to a nearest-neighbor Tm setpoint.

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
    :type overhangs: tuple

    '''

    def __init__(self, dna_object, tm=72, min_len=10, tm_undershoot=1,
                 tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                 overhangs=''):
        # TODO: type checking
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
        self.overhangs = overhangs

    def run(self):
        template = str(self.template)
        primers = _design_primer_gene(template, self.tm, self.min_len,
                                      self.tm_undershoot, self.tm_overshoot,
                                      self.end_gc, self.tm_method,
                                      self.overhangs)
        primers_dna = [sequence.DNA(seq[0], stranded='ss') for seq in primers]
        return primers_dna


def _design_primer(seq, tm=72, min_len=10, tm_undershoot=1, tm_overshoot=3,
                   end_gc=False, tm_method='finnzymes', overhang=''):
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
    :param overhang: Append the primer to this overhang sequence.
    :type overhang: str

    '''

    # Check Tm of input sequence to see if it's already too low
    seq_tm = analysis.Tm(seq, method=tm_method).run()
    if seq_tm < tm - tm_undershoot:
        err = 'Input sequence Tm is lower than primer Tm setting'
        raise Exception(err)

    max_tm = tm + tm_overshoot
    bases = min_len
    # Trim down max length to increase efficiency
    # Pretty much impossible for an annealing sequence to need more than 90bp
    seq = seq[0:90]

    # First, generate a range of primers and Tms:
    #     range: from min_len to 'tm' + tm_overshoot
    primers = []
    tms = []
    primer_tm = 0
    primer_len = 0
    while primer_tm <= max_tm and primer_len <= 80 and primer_len != len(seq):
        new_primer = seq[0:bases]
        new_tm = analysis.Tm(new_primer, method=tm_method).run()
        primers.append(new_primer)
        tms.append(new_tm)
        bases += 1
        primer_tm = tms[-1]
        primer_len = len(new_primer)

    # Trim primer list based on tm_undershoot and end_gc
    terr = tm - tm_undershoot
    #print terr
    #print tms
    primers = [primers[i] for i, x in enumerate(tms) if x >= terr]
    tms = [x for i, x in enumerate(tms) if x >= terr]

    if end_gc:
        primers = [x for x in primers if x.endswith(('C', 'G'))]
        tms = [tms[i] for i, x in enumerate(primers) if x.endswith(('C', 'G'))]

    if not primers:
        raise Exception('No primers could be generated using these settings')

    # Find the primer closest to the set Tm
    tm_diffs = [abs(x - tm) for x in tms]
    best_index = tm_diffs.index(min(tm_diffs))
    best_primer = (primers[best_index], tms[best_index])

    if overhang:
        check_alphabet(overhang)
    best_primer = (overhang + best_primer[0], best_primer[1])

    return best_primer


def _design_primer_gene(seq, tm=72, min_len=10, tm_undershoot=1,
                        tm_overshoot=3, end_gc=False, tm_method='finnzymes',
                        overhangs=''):
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
        overhangs = ['', '']

    templates = [seq, reverse_complement(seq)]
    primer_list = []

    for i, sequence in enumerate(templates):
        primer = _design_primer(sequence, tm=tm, min_len=min_len,
                                tm_undershoot=tm_undershoot,
                                tm_overshoot=tm_overshoot, end_gc=end_gc,
                                tm_method=tm_method,
                                overhang=overhangs[i])
        primer_list.append(primer)
    return primer_list
