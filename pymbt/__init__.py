'''python molecular biology tool, a cloning-oriented DNA design module.'''
try:
    import matplotlib
except ImportError:
    print "Failed to import matplotlib. Plotting sequencing won't work."

from . import analysis
from . import constants
from . import database
from . import design
from . import reaction
from . import seqio
from . import simulation
from ._sequence._dna import DNA
from ._sequence._rna import RNA
from ._sequence._peptide import Peptide
from ._sequence._dna import Primer
from ._sequence._dna import RestrictionSite
from ._sequence._dna import Feature
