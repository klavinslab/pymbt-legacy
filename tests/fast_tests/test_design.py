'''Tests primer design module.'''

from nose.tools import assert_equals
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
    primer, tm = design.DesignPrimer(dna_seq, tm=72, min_len=10,
                                     tm_undershoot=1, tm_overshoot=3,
                                     end_gc=False, tm_method='finnzymes',
                                     overhang='').run()
    assert_equals(str(primer), 'atggtgagcaagggcgaggag')


def test_design_primer_gene():
    '''Test design_primer_gene function.'''

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
    primers_list = design.DesignPrimerGene(dna_seq, tm=72, min_len=10,
                                           tm_undershoot=1, tm_overshoot=3,
                                           end_gc=False, tm_method='finnzymes',
                                           overhangs='').run()
    primers = [str(x[0]) for x in primers_list]
    assert_equals(primers, ['atggtgagcaagggcgaggag',
                            'ttacttgtacagctcgtccatgccg'])


def test_oligo_assembly():
    '''Tests output of OligoAssembly class.'''

    # Expected outputs
    olig1 = 'ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTG' + \
            'ATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAA'
    olig2 = 'TGGCCATGGAACAGGTAGTTTTCCAGTAGTGCAAATAAATTTAAGGGTAAGTTTTCCGTAT' + \
            'GTTGCATCACCTTCACCCTCTCCACTGACAGAAAATTTGTG'
    olig3 = 'TGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGC' + \
            'TTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAA'
    olig4 = 'CGTGTCTTGTAGTTCCCGTCATCTTTGAAAAATATAGTTCTTTCCTGTACATAACCTTCGG' + \
            'GCATGGCACTCTTGAAAAAGTCATGCTGTTTCATATGATCTGGG'
    olig5 = 'TTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTG' + \
            'TTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACA'
    olig6 = 'TTTGTCTGCCATGATGTATACATTGTGTGAGTTATAGTTGTATTCCAATTTGTGTCCAAGA' + \
            'ATGTTTCCATCTTCTTTAAAATCAATACCTTTTAACTCGATTCTATT'
    olig7 = 'AACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTA' + \
            'ACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCA'
    olig8 = 'TTGTGTGGACAGGTAATGGTTGTCTGGTAAAAGGACAGGGCCATCGCCAATTGGAGTATTT' + \
            'TGTTGATAATGGTCTGCTAGTTGAACGCTTCCATCTTCAATGT'
    olig9 = 'CCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAG' + \
            'ACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGA'
    olig10 = 'TTAAGCTACTAAAGCGTAGTTTTCGTCGTTTGCAGCAGGCCTTTTGTATAGTTCATCCAT' + \
             'GCCATGTGTAATCCCAGCAGCTGTTACAAACTCAAGAAGG'

    reference_oligos = [olig1, olig2, olig3, olig4, olig5, olig6, olig7, olig8,
                        olig9, olig10]
    reference_tms = [73.513413945987, 72.73367624289529, 73.73563193690501,
                     72.70706564878304, 72.72193323127516, 72.23050918438167,
                     72.07546311550084, 72.27046461560099, 73.67230272019742]

    # Run oligo synthesis on BBa_K082003
    seq = 'atgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgat' + \
          'gttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaactt' + \
          'acccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactact' + \
          'ttcggttatggtgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttc' + \
          'aagagtgccatgcccgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaac' + \
          'tacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaa' + \
          'ggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactca' + \
          'cacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattaga' + \
          'cacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggc' + \
          'gatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagat' + \
          'cccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacat' + \
          'ggcatggatgaactatacaaaaggcctgctgcaaacgacgaaaactacgctttagtagcttaa'
    dna_seq = sequence.DNA(seq)
    assembly = design.OligoAssembly(dna_seq,
                                    tm=72,
                                    length_range=(120, 120),
                                    require_even=True,
                                    start_5=True)
    assembly.run()

    # Prepare outputs vs reference
    output_oligos = [str(oligo).lower() for oligo in assembly.oligos]
    reference_oligos = [oligo.lower() for oligo in reference_oligos]

    assert_equals(output_oligos, reference_oligos)
    assert_equals(assembly.overlap_tms, reference_tms)
