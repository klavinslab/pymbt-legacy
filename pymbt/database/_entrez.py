import os
from tempfile import mkdtemp
from Bio import Entrez
import pymbt.sequence


# FIXME: If a Genome is a data structure, it should be in the DNA sequence
# module. Then rename the genome acquirer after Entrez / NCBI, etc.
# TODO: docstring
# TODO: Figure out why reading in the DNA is so slow and if it can be sped up
# - MG1655 takes 30-60 seconds to process into memory and pbt.DNA.
class Genome(pymbt.sequence.DNA):
    """Acquire a genome from Entrez

    """
    def __init__(self, genome_id):
        # Using a dummy email for now - does this violate NCBI guidelines?
        email = "loremipsum@gmail.com"
        Entrez.email = email

        print "Downloading Genome..."
        handle = Entrez.efetch(db="nucleotide", id=str(genome_id), rettype="gb",
                              retmode="text")
        print "Genome Downloaded..."
        tmpfile = os.path.join(mkdtemp(), "tmp.gb")
        with open(tmpfile, "w") as f:
            f.write(handle.read())
        genome = pymbt.seqio.read_dna(tmpfile)

        self._sequence = genome.top()
        self._bottom = genome.bottom()
        self.features = genome.features
        self.topology = genome.topology
        self.stranded = "ds"
