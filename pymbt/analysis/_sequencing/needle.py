'''Needleman-Wunsch alignment using emboss needle.'''

from tempfile import mkdtemp
from shutil import rmtree
from Bio.Emboss.Applications import NeedleallCommandline
from Bio import AlignIO

# Important note: does not produce a 'consensus' reference sequence


def needle(reference, targets, gapopen=10, gapextend=0.5):
    '''
    Do Needleman-Wunsch alignment.

    :param reference: Reference sequence.
    :type reference: str
    :param targets: Sequence(s) to align against the reference.
    :type targets: str
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float

    '''

    # Make temporary dir
    workdir = mkdtemp()

    # If input isn't a list, make it one
    if type(targets) != list:
        targets = [targets]

    # Write input files to temp dir (fasta format)
    with open(workdir + '/ref.fasta', 'w') as ref_handle:
        ref_handle.write('>ref\n')
        ref_handle.write(reference)
    with open(workdir + '/targets.fasta', 'w') as targets_handle:
        for i, target in enumerate(targets):
            targets_handle.write('>target{}\n'.format(i + 1))
            targets_handle.write('{}\n'.format(target))

    # Set up Emboss 'needle' command
    cline = NeedleallCommandline(cmd='needleall')
    cline.gapopen = gapopen
    cline.gapextend = gapextend
    cline.bsequence = workdir + '/ref.fasta'
    cline.asequence = workdir + '/targets.fasta'
    cline.aformat = 'srspair'
    cline.outfile = workdir + '/alignments.txt'

    # Run 'needle'
    cline()

    # Grab the alignments using AlignIO
    alignments = [x for x in AlignIO.parse(workdir + '/alignments.txt',
                  'emboss')]

    # Process them into a list of tuples of form: [(ref, res1), (ref, res2)],
    # etc
    aligned_refs = [x[1].seq.tostring().upper() for x in alignments]
    aligned_res = [x[0].seq.tostring().upper() for x in alignments]
    alignments = zip(aligned_refs, aligned_res)

    # Manually grab the score (AlignIO doesn't get it for some reason)
    scores = []
    with open(workdir + '/alignments.txt', 'r') as align_file:
        for line in align_file.readlines():
            if line.startswith('# Score'):
                score = float(line.strip().lstrip('# Score: '))
                scores.append(score)

    # Delete temporary dir
    rmtree(workdir)

    return {'alignments': alignments, 'scores': scores}
