'''
Tests GeneSplitter design class.

'''

from nose.tools import assert_equals
from pymbt import design
from pymbt.sequence import dna


def test_genesplitter():
    '''Test GeneSplitter class methods.'''

    seq = ''.join(['ggggacaagtttgtacaaaaaagcaggcttcaaaatgatg',
                   'tctgcttctagattggctggtactttgattccagctatgg',
                   'ctttcttgtcttgtgttagaccagaatcttgggaaccatg',
                   'tgttgaagttgttccaaatattacttatcaatgtatggaa',
                   'ttgaatttctataaaattccagataatttgccattctcta',
                   'ctaaaaatttggatttgtctttcaatccattgagacattt',
                   'gggttcttattctttcttctctttcccagaattgcaagtt',
                   'ttggatttgtctagatgtgaaattcaaactattgaagatg',
                   'gtgcttatcaatctttgtctcatttgtctactttgatttt',
                   'gactggtaatccaattcaatctttggctttgggtgctttc',
                   'tctggtttgtcttctttgcaaaaattggttgctgttgaaa',
                   'ctaatttggcttctttggaaaatttcccaattggtcattt',
                   'gaaaactttgaaagaattgaatgttgctcataatttgatt',
                   'caatctttcaaattgccagaatatttctctaatttgacta',
                   'atttggaacatttggatttgtcttctaataaaattcaatc',
                   'tatttattgtactgatttgagagttttgcatcaaatgcca',
                   'ttgttgaatttgtctttggatttgtctttgaatccaatga',
                   'atttcattcaaccaggtgctttcaaagaaattagattgca',
                   'taaattgactttgagaaataatttcgattctttgaatgtt',
                   'atgaaaacttgtattcaaggtttggctggtttggaagttc',
                   'atagattggttttgggtgaattcagaaatgaaggtaattt',
                   'ggaaaaattcgataaatctgctttggaaggtttgtgtaat',
                   'ttgactattgaagaattcagattggcttatttggattatt',
                   'atttggatgatattattgatttgttcaattgtttgactaa',
                   'tgtttcttctttctctttggtttctgttactattgaaaga',
                   'gttaaagatttctcttataatttcggttggcaacatttgg',
                   'aattggttaattgtaaattcggtcaattcccaactttgaa',
                   'attgaaatctttgaaaagattgactttcacttctaataaa',
                   'ggtggtaatgctttctctgaagttgatttgccatctttgg',
                   'aattcttggatttgtctagaaatggtttgtctttcaaagg',
                   'ttgttgttctcaatctgatttcggtactacttctttgaaa',
                   'tatttggatttgtctttcaatggtgttattactatgtctt',
                   'ctaatttcttgggtttggaacaattggaacatttggattt',
                   'ccaacattctaatttgaaacaaatgtctgaattctctgtt',
                   'ttcttgtctttgagaaatttgatttatttggatatttctc',
                   'atactcatactagagttgctttcaatggtattttcaatgg',
                   'tttgtcttctttggaagttttgaaaatggctggtaattct',
                   'ttccaagaaaatttcttgccagatattttcactgaattga',
                   'gaaatttgactttcttggatttgtctcaatgtcaattgga',
                   'acaattgtctccaactgctttcaattctttgtcttctttg',
                   'caagttttgaatatgtctcataataatttcttctctttgg',
                   'atactttcccatataaatgtttgaattctttgcaagtttt',
                   'ggattattctttgaatcatattatgacttctaaaaaacaa',
                   'gaattgcaacatttcccatcttctttggctttcttgaatt',
                   'tgactcaaaatgatttcgcttgtacttgtgaacatcaatc',
                   'tttcttgcaatggattaaagatcaaagacaattgttggtt',
                   'gaagttgaaagaatggaatgtgctactccatctgataaac',
                   'aaggtatgccagttttgtctttgaatattacttgtcaaat',
                   'gaataaaactattattggtgtttctgttttgtctgttttg',
                   'gttgtttctgttgttgctgttttggtttataaattctatt',
                   'tccatttgatgttgttggctggttgtattaaatatggtag',
                   'aggtgaaaatatttatgatgctttcgttatttattcttct',
                   'caagatgaagattgggttagaaatgaattggttaaaaatt',
                   'tggaagaaggtgttccaccattccaattgtgtttgcatta',
                   'tagagatttcattccaggtgttgctattgctgctaatatt',
                   'attcatgaaggtttccataaatctagaaaagttattgttg',
                   'ttgtttctcaacatttcattcaatctagatggtgtatttt',
                   'cgaatatgaaattgctcaaacttggcaattcttgtcttct',
                   'agagctggtattattttcattgttttgcaaaaagttgaaa',
                   'aaactttgttgagacaacaagttgaattgtatagattgtt',
                   'gtctagaaatacttatttggaatgggaagattctgttttg',
                   'ggtagacatattttctggagaagattgagaaaagctttgt',
                   'tggatggtaaatcttggaatccagaaggtactgttggtac',
                   'tggttgtaattggcaagaagctacttctatttaacaccca',
                   'gctttcttgtacaaagtggtcccc'])
    dna_seq = dna.DNA(seq)
    gs_instance = design.GeneSplitter(dna_seq)
    split = gs_instance.split(max_len=1100, min_context=200, core=60,
                              context=90, step=10, force_exhaustive=False)

    assert_equals(split['overlaps'], [(903, 963), (1933, 1993)])
    assert_equals(split['scores'], [0.7900144166666668, 0.9236661666666666])

    seq1 = 'GGGGACAAGTTTGTACAAAAAAGCAGGCTTCAAAATGATGTCTGCTTCTAGATTGGCTGGTA' + \
           'CTTTGATTCCAGCTATGGCTTTCTTGTCTTGTGTTAGACCAGAATCTTGGGAACCATGTGTT' + \
           'GAAGTTGTTCCAAATATTACTTATCAATGTATGGAATTGAATTTCTATAAAATTCCAGATAA' + \
           'TTTGCCATTCTCTACTAAAAATTTGGATTTGTCTTTCAATCCATTGAGACATTTGGGTTCTT' + \
           'ATTCTTTCTTCTCTTTCCCAGAATTGCAAGTTTTGGATTTGTCTAGATGTGAAATTCAAACT' + \
           'ATTGAAGATGGTGCTTATCAATCTTTGTCTCATTTGTCTACTTTGATTTTGACTGGTAATCC' + \
           'AATTCAATCTTTGGCTTTGGGTGCTTTCTCTGGTTTGTCTTCTTTGCAAAAATTGGTTGCTG' + \
           'TTGAAACTAATTTGGCTTCTTTGGAAAATTTCCCAATTGGTCATTTGAAAACTTTGAAAGAA' + \
           'TTGAATGTTGCTCATAATTTGATTCAATCTTTCAAATTGCCAGAATATTTCTCTAATTTGAC' + \
           'TAATTTGGAACATTTGGATTTGTCTTCTAATAAAATTCAATCTATTTATTGTACTGATTTGA' + \
           'GAGTTTTGCATCAAATGCCATTGTTGAATTTGTCTTTGGATTTGTCTTTGAATCCAATGAAT' + \
           'TTCATTCAACCAGGTGCTTTCAAAGAAATTAGATTGCATAAATTGACTTTGAGAAATAATTT' + \
           'CGATTCTTTGAATGTTATGAAAACTTGTATTCAAGGTTTGGCTGGTTTGGAAGTTCATAGAT' + \
           'TGGTTTTGGGTGAATTCAGAAATGAAGGTAATTTGGAAAAATTCGATAAATCTGCTTTGGAA' + \
           'GGTTTGTGTAATTTGACTATTGAAGAATTCAGATTGGCTTATTTGGATTATTATTTGGATGA' + \
           'TATTATTGATTTGTTCAATTGTTTGACTAATGT'
    seq2 = 'GGCTTATTTGGATTATTATTTGGATGATATTATTGATTTGTTCAATTGTTTGACTAATGTTT' + \
           'CTTCTTTCTCTTTGGTTTCTGTTACTATTGAAAGAGTTAAAGATTTCTCTTATAATTTCGGT' + \
           'TGGCAACATTTGGAATTGGTTAATTGTAAATTCGGTCAATTCCCAACTTTGAAATTGAAATC' + \
           'TTTGAAAAGATTGACTTTCACTTCTAATAAAGGTGGTAATGCTTTCTCTGAAGTTGATTTGC' + \
           'CATCTTTGGAATTCTTGGATTTGTCTAGAAATGGTTTGTCTTTCAAAGGTTGTTGTTCTCAA' + \
           'TCTGATTTCGGTACTACTTCTTTGAAATATTTGGATTTGTCTTTCAATGGTGTTATTACTAT' + \
           'GTCTTCTAATTTCTTGGGTTTGGAACAATTGGAACATTTGGATTTCCAACATTCTAATTTGA' + \
           'AACAAATGTCTGAATTCTCTGTTTTCTTGTCTTTGAGAAATTTGATTTATTTGGATATTTCT' + \
           'CATACTCATACTAGAGTTGCTTTCAATGGTATTTTCAATGGTTTGTCTTCTTTGGAAGTTTT' + \
           'GAAAATGGCTGGTAATTCTTTCCAAGAAAATTTCTTGCCAGATATTTTCACTGAATTGAGAA' + \
           'ATTTGACTTTCTTGGATTTGTCTCAATGTCAATTGGAACAATTGTCTCCAACTGCTTTCAAT' + \
           'TCTTTGTCTTCTTTGCAAGTTTTGAATATGTCTCATAATAATTTCTTCTCTTTGGATACTTT' + \
           'CCCATATAAATGTTTGAATTCTTTGCAAGTTTTGGATTATTCTTTGAATCATATTATGACTT' + \
           'CTAAAAAACAAGAATTGCAACATTTCCCATCTTCTTTGGCTTTCTTGAATTTGACTCAAAAT' + \
           'GATTTCGCTTGTACTTGTGAACATCAATCTTTCTTGCAATGGATTAAAGATCAAAGACAATT' + \
           'GTTGGTTGAAGTTGAAAGAATGGAATGTGCTACTCCATCTGATAAACAAGGTATGCCAGTTT' + \
           'TGTCTTTGAATATTACTTGTCAAATGAATAAAACTATTATTGGTGTTTCTGTTTTGTCTGTT' + \
           'TTGGTTGTTTCTGTTGTTGCTGTTTTGGTTTATAAA'
    seq3 = 'ATTGGTGTTTCTGTTTTGTCTGTTTTGGTTGTTTCTGTTGTTGCTGTTTTGGTTTATAAATT' + \
           'CTATTTCCATTTGATGTTGTTGGCTGGTTGTATTAAATATGGTAGAGGTGAAAATATTTATG' + \
           'ATGCTTTCGTTATTTATTCTTCTCAAGATGAAGATTGGGTTAGAAATGAATTGGTTAAAAAT' + \
           'TTGGAAGAAGGTGTTCCACCATTCCAATTGTGTTTGCATTATAGAGATTTCATTCCAGGTGT' + \
           'TGCTATTGCTGCTAATATTATTCATGAAGGTTTCCATAAATCTAGAAAAGTTATTGTTGTTG' + \
           'TTTCTCAACATTTCATTCAATCTAGATGGTGTATTTTCGAATATGAAATTGCTCAAACTTGG' + \
           'CAATTCTTGTCTTCTAGAGCTGGTATTATTTTCATTGTTTTGCAAAAAGTTGAAAAAACTTT' + \
           'GTTGAGACAACAAGTTGAATTGTATAGATTGTTGTCTAGAAATACTTATTTGGAATGGGAAG' + \
           'ATTCTGTTTTGGGTAGACATATTTTCTGGAGAAGATTGAGAAAAGCTTTGTTGGATGGTAAA' + \
           'TCTTGGAATCCAGAAGGTACTGTTGGTACTGGTTGTAATTGGCAAGAAGCTACTTCTATTTA' + \
           'ACACCCAGCTTTCTTGTACAAAGTGGTCCCC'
    # Prepare sequences to ensure comparability
    output_seqs = [str(sequence).lower() for sequence in
                   split['split_sequences']]
    reference_seqs = [seq1.lower(), seq2.lower(), seq3.lower()]
    assert_equals(output_seqs, reference_seqs)
