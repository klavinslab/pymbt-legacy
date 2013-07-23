'''Common Restriction Sites'''
from pymbt.sequence import DNA
from pymbt.sequence import RestrictionSite


AflII = RestrictionSite(DNA('CTTAAG'), (1, 5), name='AflII')
AgeI = RestrictionSite(DNA('ACCGGT'), (1, 5), name='AgeI')
BamHI = RestrictionSite(DNA('GGATCC'), (1, 5), name='BamHI')
DpnI = RestrictionSite(DNA('GATC'), (2, 2), name='DpnI')
DraI = RestrictionSite(DNA('TTTAAA'), (3, 3), name='DraI')
EcoRI = RestrictionSite(DNA('GAATTC'), (1, 5), name='EcoRI')
EcoRV = RestrictionSite(DNA('GATATC'), (3, 3), name='EcoRV')
FokI = RestrictionSite(DNA('GGATG'), (14, 18), name='FokI')
HindIII = RestrictionSite(DNA('AAGCTT'), (1, 5), name='HindIII')
NcoI = RestrictionSite(DNA('CCATGG'), (1, 5), name='NcoI')
NheI = RestrictionSite(DNA('GCTAGC'), (1, 5), name='NheI')
NruI = RestrictionSite(DNA('TCGCGA'), (3, 3), name='NruI')
PmeI = RestrictionSite(DNA('GTTTAAAC'), (4, 4), name='PmeI')
PstI = RestrictionSite(DNA('CTGCAG'), (5, 1), name='PstI')
SpeI = RestrictionSite(DNA('ACTAGT'), (1, 5), name='SpeI')
XbaI = RestrictionSite(DNA('TCTAGA'), (1, 5), name='XbaI')
XhoI = RestrictionSite(DNA('CTCGAG'), (1, 5), name='XhoI')
XmaI = RestrictionSite(DNA('CCCGGG'), (1, 5), name='XmaI')
SnaBI = RestrictionSite(DNA('TACGTA'), (3, 3), name='SnaBI')
AclI = RestrictionSite(DNA('AACGTT'), (2, 3), name='AclI')
