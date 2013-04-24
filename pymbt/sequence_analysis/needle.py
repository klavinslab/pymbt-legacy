'''Needleman-Wunsch alignment using emboss needle.'''

from tempfile import mkdtemp
from shutil import rmtree
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Emboss.Applications import NeedleallCommandline
from Bio import AlignIO


def needle(seq1, seq2):
    '''
    Do Needleman-Wunsch alignment.

    :param seq1: First sequence.
    :type seq1: str
    :param seq2: Second sequence.
    :type seq2: str

    '''

    workdir = mkdtemp()
    aseq_handle = open(workdir + '/aseq.fasta', 'w')
    aseq_handle.write('>seq1\n')
    aseq_handle.write(seq1)
    aseq_handle.close()
    bseq_handle = open(workdir + '/bseq.fasta', 'w')
    bseq_handle.write('>seq2\n')
    bseq_handle.write(seq2)
    bseq_handle.close()

    cline = NeedleCommandline(cmd='needle')
    cline.gapopen = 10
    cline.gapextend = 0.5
    cline.asequence = workdir + '/aseq.fasta'
    cline.bsequence = workdir + '/bseq.fasta'
    cline.outfile = workdir + '/needle.txt'
    cline()
    align = AlignIO.read(workdir + '/needle.txt', 'emboss')
    align_file = open(workdir + '/needle.txt', 'r')
    handle = align_file.readlines()
    score = [x for x in handle if x[0:7] == '# Score']
    score = float(score[0].strip()[9:])

    rmtree(workdir)
    return align, score


def needleall(seq1, seq2s):
    '''
    Do Needleman-Wunsch alignment using EMBOSS NeedleAll.

    :param seq1: First sequence.
    :type seq1: str
    :param seq2s: List of second sequences to align with the seq1.
    :type seq2s: list

    '''

    workdir = mkdtemp()
    aseq_handle = open(workdir + '/aseq.fasta', 'w')
    aseq_handle.write('>seq1\n')
    aseq_handle.write(seq1)
    aseq_handle.close()
    bseq_handle = open(workdir + '/bseq.fasta', 'w')
    if type(seq2s) == str:
        seq2s = [seq2s]
    elif type(seq2s) == list:
        pass
    else:
        raise ValueError('second input must be list of sequences')
    for i, sequence in enumerate(seq2s):
        bseq_handle.write('>seq' + str(i + 2) + '\n')
        bseq_handle.write(sequence + '\n')
    bseq_handle.close()

    cline = NeedleallCommandline(cmd='needleall')
    cline.gapopen = 10
    cline.gapextend = 0.5
    cline.asequence = workdir + '/aseq.fasta'
    cline.bsequence = workdir + '/bseq.fasta'
    cline.outfile = workdir + '/needle.txt'
    cline.aformat = 'srspair'
    cline()
    align_handle = AlignIO.parse(workdir + '/needle.txt', 'emboss')
    align = [x for x in align_handle]
    align_file = open(workdir + '/needle.txt', 'r')
    handle = align_file.readlines()
    score = [x for x in handle if x[0:7] == '# Score']
    if len(score) > 1:
        score = [float(x.strip()[9:]) for x in score]
    else:
        score = float(score[0].strip()[9:])
    to_return = [(align[i], score[i]) for i, x in enumerate(align)]

    rmtree(workdir)
    return to_return
