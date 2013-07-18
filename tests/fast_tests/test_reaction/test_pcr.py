'''Test functionality of PCR class of reaction module.'''

from pymbt import design, reaction, sequence
from nose.tools import assert_equal


def test_basic():
    to_amplify = 'atgtctaaaggtgaagaattattcactggtgttgtcccaatgctgctggtattacc' + \
                 'catggtattgatgaattgtacaaatag'
    template = sequence.DNA(to_amplify)
    forward, reverse = design.design_primer_pcr(template)

    amplicon = reaction.pcr(template, forward, reverse)
    assert_equal(amplicon, template)
