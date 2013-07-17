from pymbt import seqio, reaction


a = seqio.read_dna('/home/nick/gibson_test/a.ape')
a1 = seqio.read_dna('/home/nick/gibson_test/a1.fasta')
a2 = seqio.read_dna('/home/nick/gibson_test/a2.fasta')
a3 = seqio.read_dna('/home/nick/gibson_test/a3.fasta')

final = reaction.Gibson([a1, a2, a3]).run()
print 'expected plasmid size: {}'.format(len(a))
print 'fragment lengths: {}, {}, {}'.format(len(a1), len(a2), len(a3))
print 'Gibson simulation product size: {}'.format(len(final))
