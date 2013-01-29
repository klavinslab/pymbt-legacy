#Takes in two strings from an alignment and makes a report. This is just a test of writing my own mismatch/deletion/insertion report
# Note: reports deviations from first sequence (seq_reference)

def align_analyze(seq_reference,seq_result):
	# Figure out the right way to handle this error in python
	if len(seq_reference)!=len(seq_result):
		print 'Sequences should be of identical length - this report is wrong'

	# Current method: probably the slowest possible, but works. Evaluates every position as either a mismatch, deletion, or insertion

	#convert to upper and lists
	seq_reference = seq_reference.upper()
	seq_result = seq_result.upper()
	ref = list(seq_reference)
	res = list(seq_result)

	# Set up lists (indexes) of mismatches, deletions, insertions
	mismatches = []
	deletions = []
	insertions = []
	for i,v in enumerate(res):
		if v == '-':
			deletions.append(i)
		elif ref[i]=='-':
				insertions.append(i)
		elif v!=ref[i]:
				mismatches.append(i)
	# Remove initial and final string of 'deletions' - they are just beginning and end gaps
	del_diffs = []
	for i in range(len(deletions)-1):
		del_diffs.append(deletions[i+1]-deletions[i])
	gap_front = [i for i,v in enumerate(del_diffs) if v != 1][0]
	gap_rear = len(del_diffs)-[i for i,v in enumerate(reversed(del_diffs)) if v != 1][0]
	ref_nogap = ref[gap_front:gap_rear]
	res_nogap = res[gap_front:gap_rear]

	print(''.join(ref_nogap))
	print(''.join(res_nogap))
	
	#return(ref_nogap)
	#return((first_notn,last_notn))

	#print 'mismatches: ' + ''.join(["%s," % el for el in mismatches])
	#print 'deletions: ' + ''.join(["%s," % el for el in deletions])
	#print 'insertions: ' + ''.join(["%s," % el for el in insertions])
