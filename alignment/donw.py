# Does needleman-wunsch alignment using emboss needle writing to temp file, returns alignment object

from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

def donw(seq1, seq2):
	cline = NeedleCommandline(cmd='needle')
	cline.gapopen = 10
	cline.gapextend = 0.5
	cline.asequence = 'asis:'+seq1
	cline.bsequence = 'asis:'+seq2
	cline.outfile = '/tmp/needle.txt'
	cline()
	align = AlignIO.read('/tmp/needle.txt','emboss')
	return(align)
