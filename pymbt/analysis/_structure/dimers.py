'''Check for primer dimers using Nupack.'''
from pymbt.analysis import Nupack


def dimers(primer1, primer2, concentrations=[5e-7, 5e-7, 3e-11]):
    '''Check for probability of primer dimers given primer and template concs.

    :param primer1: Forward primer
    :type primer1: pymbt.sequence.DNA
    :param primer2: Reverse primer
    :type primer2: pymbt.sequence.DNA
    :param template: DNA template
    :type template: pymbt.sequence.DNA
    :param concentrations: list of concentrations for primer1, primer2,
                           template. Defaults are those for PCR with 1kb
                           template.
    :type concentrations: list

    '''
    # It is not reasonable (yet) to use a long template for doing these
    # computations directly, as NUPACK does an exhaustive search (and does not
    # make use of multiple cores) and takes far too long. Instead, this
    # function just compares self-self vs. self-other binding

    # Simulate binding of template vs. primers
    nupack = Nupack([primer1.primer(), primer2.primer,
                     primer1.primer().reverse_complement(),
                     primer2.primer().reverse_complement()])
    # Include reverse complement concentration
    concentrations.append(concentrations[2])
    nupack_concs = nupack.concentrations(2, conc=concentrations)
    dimer_conc = nupack_concs['concentrations'][5]
    primer1_template = nupack_concs['concentrations'][6]
    primer2_template = nupack_concs['concentrations'][10]
    return primer1_template, primer2_template, dimer_conc
