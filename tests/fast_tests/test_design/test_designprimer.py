'''Tests primer design module.'''
import warnings
from nose.tools import assert_equals, assert_not_equal, assert_raises
from pymbt import design
from pymbt import sequence


def test_design_primer():
    '''Test design_primer function.'''
    seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggc' + \
          'gacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaag' + \
          'ctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgacc' + \
          'accttcggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttc' + \
          'ttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggc' + \
          'aactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctg' + \
          'aagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaac' + \
          'agccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatc' + \
          'cgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatc' + \
          'ggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaa' + \
          'gaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcact' + \
          'ctcggcatggacgagctgtacaagtaa'
    dna_seq = sequence.DNA(seq)
    primer = design.design_primer(dna_seq, tm=72, min_len=10,
                                  tm_undershoot=1, tm_overshoot=3,
                                  end_gc=False, tm_parameters='cloning',
                                  overhang=None)
    assert_equals(str(primer), 'atggtgagcaagggcgaggag')
    # Ensure that overhang is appropriately applied
    overhang_primer = design.design_primer(dna_seq, tm=72, min_len=10,
                                           tm_undershoot=1, tm_overshoot=3,
                                           end_gc=False,
                                           tm_parameters='cloning',
                                           overhang=sequence.DNA('GATCGATAT'))
    assert_equals(str(overhang_primer), 'gatcgatatatggtgagcaagggcgaggag')
    # If sequence is too short (too low of Tm), raise ValueError
    too_short = sequence.DNA('at')
    assert_raises(ValueError, design.design_primer, too_short, tm=72)
    # Should design different primers (sometimes) if ending on GC is preferred
    diff_template = sequence.DNA('GATCGATCGATACGATCGATATGCGATATGATCGATAT')
    nogc = design.design_primer(diff_template, tm=72, min_len=10,
                                tm_undershoot=1, tm_overshoot=3, end_gc=False,
                                tm_parameters='cloning', overhang=None)
    withgc = design.design_primer(diff_template, tm=72, min_len=10,
                                  tm_undershoot=1, tm_overshoot=3,
                                  end_gc=True, tm_parameters='cloning',
                                  overhang=None)
    assert_not_equal(nogc, withgc)
    # Should raise ValueError if it's impossible to create an end_gc primer
    end_at_template = sequence.DNA('ATGCGATACGATACGCGATATGATATATatatatat' +
                                   'ATAAaaaaaaaaaattttttttTTTTTTTTTTTTTT' +
                                   'TTTTTTTTTT')
    assert_raises(ValueError, design.design_primer, end_at_template,
                  end_gc=True, tm=72)
    # If there's structure, should issue a warning
    structure_template = sequence.DNA('ATGCGATCGATAGGCGA')
    structure_template += structure_template.reverse_complement()
    with warnings.catch_warnings(True) as w:
        design.design_primer(structure_template, structure=True, tm=72)
        assert len(w) > 0


def test_design_primers():
    '''Test design_primers function.'''
    seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggc' + \
          'gacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaag' + \
          'ctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgacc' + \
          'accttcggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttc' + \
          'ttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggc' + \
          'aactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctg' + \
          'aagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaac' + \
          'agccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatc' + \
          'cgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatc' + \
          'ggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaa' + \
          'gaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcact' + \
          'ctcggcatggacgagctgtacaagtaa'
    dna_seq = sequence.DNA(seq)
    primers_list = design.design_primers(dna_seq, tm=72, min_len=10,
                                         tm_undershoot=1, tm_overshoot=3,
                                         end_gc=False,
                                         tm_parameters='cloning',
                                         overhangs=None)
    primers = [str(x.primer()) for x in primers_list]
    assert_equals(primers, ['atggtgagcaagggcgaggag',
                            'ttacttgtacagctcgtccatgccg'])
