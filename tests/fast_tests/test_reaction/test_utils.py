'''
Tests utils submodule of reaction module.

'''

from nose.tools import assert_equal, assert_raises
from pymbt import reaction
from pymbt import sequence


def test_convert_sequence():
    '''Tests DNA translation function.'''

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
    dna = sequence.DNA(seq)
    prot = 'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLV' + \
           'TTFGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRI' + \
           'ELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQN' + \
           'TPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK'
    prot = prot.lower()
    rna = reaction.utils.convert_sequence(dna, 'rna')
    r_trans = reaction.utils.convert_sequence(rna, 'dna')
    trans = reaction.utils.convert_sequence(rna, 'peptide')
    assert_equal(str(trans), prot)
    assert_equal(str(r_trans), seq)

    # Gapped sequence shouldfail
    def convert_gap_seq(seq):
        reaction.utils.convert_sequence(seq, 'rna')

    gapped_dna = sequence.DNA('atg-')

    assert_raises(ValueError, convert_gap_seq, gapped_dna)

    # Sequence without stop codon should still work
    nostop_dna = sequence.DNA('atgaaaaaaaaaaaa')
    nostop_rna = reaction.utils.convert_sequence(nostop_dna, 'rna')
    nostop_peptide = reaction.utils.convert_sequence(nostop_rna, 'peptide')
    assert_equal(str(nostop_rna), 'augaaaaaaaaaaaa')
    assert_equal(str(nostop_peptide), 'mkkkk')
