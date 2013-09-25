'''Needleman-Wunsch alignment using emboss needle.'''
import os
from tempfile import mkdtemp
from shutil import rmtree
from pymbt import sequence
from Bio.Emboss.Applications import NeedleallCommandline
from Bio import AlignIO
from Bio.Application import ApplicationError


def needle(reference, target, gapopen=10, gapextend=0.5):
    '''Do Needleman-Wunsch alignment.

    :param reference: Reference sequence.
    :type reference: pymbt.sequence.DNA
    :param target: Sequence to align against the reference.
    :type target: pymbt.sequence.DNA
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float

    '''
    # Check to see if 'needleall' command is installed - if not, give useful
    # warning
    msg = "The command 'needleall' is not available. You may need to " + \
          "install EMBOSS."
    needleall_executables = []
    for path in os.environ['PATH'].split(os.pathsep):
        exepath = os.path.join(path, 'needleall')
        needleall_executables.append(os.access(exepath, os.X_OK))
    if not any(needleall_executables):
        raise Exception(msg)

    # Make temporary dir and move there
    workdir = mkdtemp()
    old_dir = os.getcwd()
    # Purpose of moving is so that 'needleall.error' doesn't pollute other dirs
    os.chdir(workdir)

    # Write input files to temp dir (fasta format)
    with open(workdir + '/reference.fasta', 'w') as reference_handle:
        reference_handle.write('>reference\n')
        reference_handle.write('{}\n'.format(str(reference)))
    with open(workdir + '/target.fasta', 'w') as target_handle:
        target_handle.write('>target\n')
        target_handle.write('{}\n'.format(str(target)))

    # Set up Emboss 'needle' command
    cline = NeedleallCommandline(cmd='needleall')
    cline.gapopen = gapopen
    cline.gapextend = gapextend
    cline.bsequence = workdir + '/reference.fasta'
    cline.asequence = workdir + '/target.fasta'
    cline.aformat = 'srspair'
    cline.outfile = workdir + '/alignment.txt'

    # Run 'needle'
    try:
        cline()
    except ApplicationError:
        # TODO: replace this with a useful message
        os.chdir(old_dir)
        rmtree(workdir)
        raise AlignmentError('Failed to align sequences.')

    # Grab the alignment using AlignIO
    alignio = AlignIO.read(workdir + '/alignment.txt', 'emboss')

    # Process them into a list of tuples of form: [(ref, res1), (ref, res2)],
    # etc
    aligned_reference = sequence.DNA(alignio[1].seq.tostring())
    aligned_result = sequence.DNA(alignio[0].seq.tostring())

    # Manually grab the score (AlignIO doesn't get it for some reason)
    with open(workdir + '/alignment.txt', 'r') as align_file:
        for line in align_file.readlines():
            if line.startswith('# Score'):
                score = float(line.strip().lstrip('# Score: '))
    if not score:
        raise Exception("Could not find alignment score (needle)!")

    # Leave and delete temporary dir
    os.chdir(old_dir)
    rmtree(workdir)

    # 'alignment' is a list. Each element of the list is a tuple of (1)
    # the aligned reference sequence and (2) the aligned target sequene
    return aligned_reference, aligned_result, score


class AlignmentError(Exception):
    '''Exception to throw when an alignment fails to run at all.'''
    pass
