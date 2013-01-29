import sys
from Bio import SeqIO
from donw import donw as donw
from pymbt.alignment.alignanalyze import align_analyze as analyze

def alignreport(path1,path2):
	ref = SeqIO.read(path1,'fasta').seq.tostring()
	res = SeqIO.read(path2,'fasta').seq.tostring()

	# Trim out Ns, keep largest non-N contiguous region
	newref = max(ref.split("N"),key=len)
	newres = max(res.split("N"),key=len)
	aligned = donw(newref,newres)
	return(aligned)
	aligned0 = aligned[0,:].seq.tostring()
	aligned1 = aligned[1,:].seq.tostring()
	print(analyze(aligned0,aligned1))
#	return(analyze(aligned0,aligned1))

if __name__=="__main__":
	if len(sys.argv) == 3:
		alignreport(sys.argv[1],sys.argv[2])
#		msg = plotabi(sys.argv[1])
#		print(msg)
	else:
		print("No input?")
