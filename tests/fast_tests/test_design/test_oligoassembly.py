'''
Tests for the OligoAssembly design class.

'''

from nose.tools import assert_equals
from pymbt import design
from pymbt import sequence


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
    reference_tms = [73.513413945987, 72.73367624289534, 73.73563193690484,
                     72.70706564878299, 72.72193323127533, 72.23050918438184,
                     72.07546311550101, 72.27046461560093, 73.67230272019759]

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

    print [abs(x / y - 1) for x, y in zip(assembly.overlap_tms, reference_tms)]
    print assembly.overlap_tms
    print reference_tms
    assert_equals(output_oligos, reference_oligos)
    assert_equals(assembly.overlap_tms, reference_tms)
