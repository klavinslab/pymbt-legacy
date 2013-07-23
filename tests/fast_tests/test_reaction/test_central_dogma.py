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
    transcription_output = reaction.transcribe(test_dna)
    assert_equal(transcription_output, reference_rna)

    # Coding RNA should exclude anything after a stop codon
    coding_rna_output = reaction.coding_sequence(transcription_output)
    assert_equal(coding_rna_output, reference_rna[:-3])

    # Should fail is sequence lacks start codon or stop codon
    assert_raises(ValueError, reaction.coding_sequence,
                  reaction.transcribe(sequence.DNA('aaatag')))
    assert_raises(ValueError, reaction.coding_sequence,
                  reaction.transcribe(sequence.DNA('atgaaa')))


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
    translation_output = reaction.translate(test_rna)
    assert_equal(translation_output, reference_peptide)

    # Coding peptide should exclude anything after a stop codon
    coding_rna = reaction.coding_sequence(test_rna)
    coding_peptide = reaction.translate(coding_rna)
    assert_equal(coding_peptide, reference_peptide)


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
    r_transcription = reaction.reverse_transcribe(test_rna)
    assert_equal(r_transcription, ref_dna)
