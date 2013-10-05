'''Needleman-Wunsch alignment using emboss needle.'''
import os
from tempfile import mkdtemp
from shutil import rmtree
import multiprocessing
from pymbt import sequence
from Bio.Emboss.Applications import NeedleallCommandline
from Bio import AlignIO
from Bio.Application import ApplicationError


class AlignmentError(Exception):
    '''Exception to throw when an alignment fails to run at all.'''
    pass


def needle(reference, query, gapopen=10, gapextend=0.5):
    '''Do Needleman-Wunsch alignment.

    :param reference: Reference sequence.
    :type reference: pymbt.sequence.DNA
    :param query: Sequence to align against the reference.
    :type query: pymbt.sequence.DNA
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float
    :returns: (aligned reference, aligned query, score)
    :rtype: tuple of two pymbt.sequence.DNA instances and a float
    :raises: Exception if EMBOSS can't be found

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
    with open(workdir + '/query.fasta', 'w') as query_handle:
        query_handle.write('>query\n')
        query_handle.write('{}\n'.format(str(query)))

    # Set up Emboss 'needle' command
    cline = NeedleallCommandline(cmd='needleall')
    cline.gapopen = gapopen
    cline.gapextend = gapextend
    cline.bsequence = workdir + '/reference.fasta'
    cline.asequence = workdir + '/query.fasta'
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

    # Leave and delete temporary dir
    os.chdir(old_dir)
    rmtree(workdir)

    # 'alignment' is a list. Each element of the list is a tuple of (1)
    # the aligned reference sequence and (2) the aligned query sequene
    return aligned_reference, aligned_result, score


def run_needle(args):
    """Run needle command using 4-tuple of the arguments (in the same order)
    as is used for needle. Necessary to make picklable function for
    multiprocessing."""
    return needle(*args)


def needle_multiprocessing(references, queries, gapopen=10, gapextend=0.5):
    """Batch process of sequencing split over several cores. Acts just like
    needle but sequence inputs are lists.

    :param references: References sequence.
    :type references: pymbt.sequence.DNA list
    :param queries: Sequences to align against the reference.
    :type queries: pymbt.sequence.DNA list
    :param gapopen: Penalty for opening a gap.
    :type gapopen: float
    :param gapextend: Penalty for extending a gap.
    :type gapextend: float
    :returns: a list of the same output as pymbt.sequence.needle
    :rtype: list

    """
    pool = multiprocessing.Pool()
    try:
        args_list = [list(x) + [gapopen, gapextend] for x in
                     zip(references, queries)]
        aligned = pool.map(run_needle, args_list)
    except KeyboardInterrupt:
        pool.terminate()
        raise KeyboardInterrupt

    return aligned
