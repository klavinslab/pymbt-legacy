from nose.tools import assert_equals
from pymbt.oligo_synthesis import oligo_assembly


def test_oligo_assembly():
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

    expected_oligos = [olig1, olig2, olig3, olig4, olig5, olig6, olig7, olig8,
                       olig9, olig10]
    expected_tms = [73.513413945987, 72.73367624289529, 73.73563193690501,
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
    designed = oligo_assembly.OligoAssembly(seq,
                                            tm=72,
                                            oligo_size=120,
                                            require_even=True,
                                            start_5=True)

    assert_equals(designed.oligos, expected_oligos)
    assert_equals(designed.overlap_tms, expected_tms)
