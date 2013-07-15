'''Test functionality of PCR class of reaction module.'''

from pymbt import design, reaction, sequence
from nose.tools import assert_equal


def test_basic():
    to_amplify = 'atgtctaaaggtgaagaattattcactggtgttgtcccaatgctgctggtattacc' + \
                 'catggtattgatgaattgtacaaatag'
    template = sequence.DNA(to_amplify)
    forward = design.DesignPrimer(template).run()
    reverse = design.DesignPrimer(template.reverse_complement()).run()

    pcr = reaction.PCR(forward, reverse, template)
    amplicon = pcr.run()
    assert_equal(amplicon, template)
