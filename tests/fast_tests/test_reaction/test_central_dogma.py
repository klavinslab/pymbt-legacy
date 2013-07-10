'''
Tests for central dogma submodule of reaction module.

'''

from nose.tools import assert_equal, assert_raises
from pymbt import sequence
from pymbt import reaction


def test_transcription():
    test_dna = sequence.DNA('ATGATGGGCAGTGTCGAATTAAATCTGCGTGAGACAGAATTGTGTT' +
                            'TGGGACTACCAGGCGGTGATACAGTTGCACCAGTAACAGGAAACAA' +
                            'AAGAGGATTCTCTGAAACAGTAGATTTGAAACTTAATTTGAACAAT' +
                            'GAGCCAGCCAACAAGGAAGGTTCCACCACTCATGACGTCGTCACAT' +
                            'TTGATAGTAAAGAAAAGAGTGCGTGTCCAAAAGATCCAGCTAAGCC' +
                            'ACCTGCCAAGGCTCAAGTCGTCGGATGGCCACCTGTGAGATCTTAT' +
                            'AGAAAGAACGTAATGGTTTCTTGTCAGAAGTCCAGTGGTGGTCCTG' +
                            'AAGCAGCGGCTtgaaaa')
    reference_rna = sequence.RNA('AUGAUGGGCAGUGUCGAAUUAAAUCUGCGUGAGACAGAAUU' +
                                 'GUGUUUGGGACUACCAGGCGGUGAUACAGUUGCACCAGUAA' +
                                 'CAGGAAACAAAAGAGGAUUCUCUGAAACAGUAGAUUUGAAA' +
                                 'CUUAAUUUGAACAAUGAGCCAGCCAACAAGGAAGGUUCCAC' +
                                 'CACUCAUGACGUCGUCACAUUUGAUAGUAAAGAAAAGAGUG' +
                                 'CGUGUCCAAAAGAUCCAGCUAAGCCACCUGCCAAGGCUCAA' +
                                 'GUCGUCGGAUGGCCACCUGUGAGAUCUUAUAGAAAGAACGU' +
                                 'AAUGGUUUCUUGUCAGAAGUCCAGUGGUGGUCCUGAAGCAG' +
                                 'CGGCUugaaaa')
    # Basic transcription should work
    transcription = reaction.Transcription(test_dna)
    transcription_output = transcription.run()
    assert_equal(transcription_output, reference_rna)

    # Coding RNA should exclude anything after a stop codon
    coding_rna_output = transcription.get_coding_rna()
    assert_equal(coding_rna_output, reference_rna[:-3])

    # If coding rna is wanted, should work even if .run() hasn't happened
    transcription_unrun = reaction.Transcription(test_dna)
    coding_rna_unrun_output = transcription_unrun.get_coding_rna()
    assert_equal(coding_rna_unrun_output, reference_rna[:-3])

    # Should fail is sequence lacks start codon or stop codon
    def nostart():
        transcription_nostart = reaction.Transcription(sequence.DNA('aaatag'))
        transcription_nostart.get_coding_rna()
    assert_raises(Exception, nostart)

    def nostop():
        transcription_nostop = reaction.Transcription(sequence.DNA('atgaaa'))
        transcription_nostop.get_coding_rna()
    assert_raises(Exception, nostop)


def test_translation():
    test_rna = sequence.RNA('AUGAUGGGCAGUGUCGAAUUAAAUCUGCGUGAGACAGAAUU' +
                            'GUGUUUGGGACUACCAGGCGGUGAUACAGUUGCACCAGUAA' +
                            'CAGGAAACAAAAGAGGAUUCUCUGAAACAGUAGAUUUGAAA' +
                            'CUUAAUUUGAACAAUGAGCCAGCCAACAAGGAAGGUUCCAC' +
                            'CACUCAUGACGUCGUCACAUUUGAUAGUAAAGAAAAGAGUG' +
                            'CGUGUCCAAAAGAUCCAGCUAAGCCACCUGCCAAGGCUCAA' +
                            'GUCGUCGGAUGGCCACCUGUGAGAUCUUAUAGAAAGAACGU' +
                            'AAUGGUUUCUUGUCAGAAGUCCAGUGGUGGUCCUGAAGCAG' +
                            'CGGCUugaaaa')
    reference_peptide = sequence.Peptide('MMGSVELNLRETELCLGLPGGDTVAPVTGNK' +
                                         'RGFSETVDLKLNLNNEPANKEGSTTHDVVTF' +
                                         'DSKEKSACPKDPAKPPAKAQVVGWPPVRSYR' +
                                         'KNVMVSCQKSSGGPEAAA')
    # Basic transcription should work
    translation = reaction.Translation(test_rna)
    translation_output = translation.run()
    assert_equal(translation_output, reference_peptide)

    # Coding peptide should exclude anything after a stop codon
    coding_peptide_output = translation.get_coding_peptide()
    assert_equal(coding_peptide_output, reference_peptide)

    # Should fail is sequence lacks start codon or stop codon
    def nostart():
        translation_nostart = reaction.Transcription(sequence.RNA('aaauag'))
        translation_nostart.get_coding_peptide()
    assert_raises(Exception, nostart)

    def nostop():
        translation_nostop = reaction.Transcription(sequence.RNA('augaaa'))
        translation_nostop.get_coding_peptide()
    assert_raises(Exception, nostop)


def test_reverse_transcription():
    test_rna = sequence.RNA('AUGAUGGGCAGUGUCGAAUUAAAUCUGCGUGAGACAGAAUU' +
                            'GUGUUUGGGACUACCAGGCGGUGAUACAGUUGCACCAGUAA' +
                            'CAGGAAACAAAAGAGGAUUCUCUGAAACAGUAGAUUUGAAA' +
                            'CUUAAUUUGAACAAUGAGCCAGCCAACAAGGAAGGUUCCAC' +
                            'CACUCAUGACGUCGUCACAUUUGAUAGUAAAGAAAAGAGUG' +
                            'CGUGUCCAAAAGAUCCAGCUAAGCCACCUGCCAAGGCUCAA' +
                            'GUCGUCGGAUGGCCACCUGUGAGAUCUUAUAGAAAGAACGU' +
                            'AAUGGUUUCUUGUCAGAAGUCCAGUGGUGGUCCUGAAGCAG' +
                            'CGGCUugaaaa')
    ref_dna = sequence.DNA('ATGATGGGCAGTGTCGAATTAAATCTGCGTGAGACAGAATTGTGTT' +
                           'TGGGACTACCAGGCGGTGATACAGTTGCACCAGTAACAGGAAACAA' +
                           'AAGAGGATTCTCTGAAACAGTAGATTTGAAACTTAATTTGAACAAT' +
                           'GAGCCAGCCAACAAGGAAGGTTCCACCACTCATGACGTCGTCACAT' +
                           'TTGATAGTAAAGAAAAGAGTGCGTGTCCAAAAGATCCAGCTAAGCC' +
                           'ACCTGCCAAGGCTCAAGTCGTCGGATGGCCACCTGTGAGATCTTAT' +
                           'AGAAAGAACGTAATGGTTTCTTGTCAGAAGTCCAGTGGTGGTCCTG' +
                           'AAGCAGCGGCTtgaaaa')

    # Basic transcription should work
    r_transcription = reaction.ReverseTranscription(test_rna)
    r_transcription_output = r_transcription.run()
    assert_equal(r_transcription_output, ref_dna)
