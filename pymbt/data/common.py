'''Defines data and parameters in an easily resuable format.'''

# Common sequence alphabets.
ALPHABETS = {
    'dna': 'ATGCNatgcn-',
    'rna': 'AUGCNaugcn',
    'pep': 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy'}

COMPLEMENTS = {
    'dna': ('ATGCNatgcn-', 'TACGNtacgn-'),
    'rna': ('AUGCNaugcn', 'UACGNuacgn')}


# The standard codon table.
CODON_TABLE = {
    'A': ['GCG', 'GCA', 'GCT', 'GCC'],
    'R': ['AGG', 'AGA', 'CGG', 'CGA', 'CGT', 'CGC'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    '*': ['TGA', 'TAG', 'TAA'],
    'Q': ['CAG', 'CAA'],
    'E': ['GAG', 'GAA'],
    'G': ['GGG', 'GGA', 'GGT', 'GGC'],
    'H': ['CAT', 'CAC'],
    'I': ['ATA', 'ATT', 'ATC'],
    'L': ['TTG', 'TTA', 'CTG', 'CTA', 'CTT', 'CTC'],
    'K': ['AAG', 'AAA'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCG', 'CCA', 'CCT', 'CCC'],
    'S': ['AGT', 'AGC', 'TCG', 'TCA', 'TCT', 'TCC'],
    'T': ['ACG', 'ACA', 'ACT', 'ACC'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTG', 'GTA', 'GTT', 'GTC']}


# Saccharomyces cerevisiae
# source: http://www.kazusa.or.jp/codon/
# (which cites GenBank, i.e. yeast genome project CDS database)
CODON_FREQ = {
    'sc': {
        'GCG': 0.109972396541529,
        'GCA': 0.288596474496094,
        'GCT': 0.377014739102356,
        'GCC': 0.224416389860021,
        'AGG': 0.208564104515562,
        'AGA': 0.481137590939125,
        'CGG': 0.0392677130215486,
        'CGA': 0.0676728924436203,
        'CGT': 0.144572019635586,
        'CGC': 0.0587856794445578,
        'AAT': 0.589705127199784,
        'AAC': 0.410294872800217,
        'GAT': 0.65037901553924,
        'GAC': 0.34962098446076,
        'TGT': 0.629812614586062,
        'TGC': 0.370187385413938,
        'TGA': 0.303094329334787,
        'TAG': 0.225736095965104,
        'TAA': 0.471169574700109,
        'CAG': 0.307418833439535,
        'CAA': 0.692581166560465,
        'GAG': 0.296739610207218,
        'GAA': 0.703260389792782,
        'GGG': 0.119057918187951,
        'GGA': 0.215422869017838,
        'GGT': 0.472217600813099,
        'GGC': 0.193301611981112,
        'CAT': 0.636710255236351,
        'CAC': 0.363289744763649,
        'ATA': 0.273331091899568,
        'ATT': 0.462925823433014,
        'ATC': 0.263743084667417,
        'TTG': 0.286319859527146,
        'TTA': 0.275534472444779,
        'CTG': 0.110440170850593,
        'CTA': 0.141277445174148,
        'CTT': 0.129115062940288,
        'CTC': 0.0573129890630467,
        'AAG': 0.423936637198697,
        'AAA': 0.576063362801303,
        'ATG': 1,
        'TTT': 0.586126603840976,
        'TTC': 0.413873396159024,
        'CCG': 0.120626895854398,
        'CCA': 0.417143753704543,
        'CCT': 0.307740315888567,
        'CCC': 0.154489034552491,
        'AGT': 0.159245398699046,
        'AGC': 0.109749229743856,
        'TCG': 0.0963590866114069,
        'TCA': 0.210157220085731,
        'TCT': 0.264456618519558,
        'TCC': 0.160032446340401,
        'ACG': 0.135583991997041,
        'ACA': 0.302413913478422,
        'ACT': 0.345237040780705,
        'ACC': 0.216765053743832,
        'TGG': 1,
        'TAT': 0.559573963633711,
        'TAC': 0.440426036366289,
        'GTG': 0.190897642582249,
        'GTA': 0.208783185960798,
        'GTT': 0.391481704636128,
        'GTC': 0.208837466820824}}

# Codon usage for S cerevisiae
CODON_FREQ_SC_NESTED = {
    'A': {'GCG': 0.109972396541529,
          'GCA': 0.288596474496094,
          'GCT': 0.377014739102356,
          'GCC': 0.224416389860021},
    'R': {'AGG': 0.208564104515562,
          'AGA': 0.481137590939125,
          'CGG': 0.0392677130215486,
          'CGA': 0.0676728924436203,
          'CGT': 0.144572019635586,
          'CGC': 0.0587856794445578},
    'N': {'AAT': 0.589705127199784,
          'AAC': 0.410294872800217},
    'D': {'GAT': 0.65037901553924,
          'GAC': 0.34962098446076},
    'C': {'TGT': 0.629812614586062,
          'TGC': 0.370187385413938},
    '*': {'TGA': 0.303094329334787,
          'TAG': 0.225736095965104,
          'TAA': 0.471169574700109},
    'Q': {'CAG': 0.307418833439535,
          'CAA': 0.692581166560465},
    'E': {'GAG': 0.296739610207218,
          'GAA': 0.703260389792782},
    'G': {'GGG': 0.119057918187951,
          'GGA': 0.215422869017838,
          'GGT': 0.472217600813099,
          'GGC': 0.193301611981112},
    'H': {'CAT': 0.636710255236351,
          'CAC': 0.363289744763649},
    'I': {'ATA': 0.273331091899568,
          'ATT': 0.462925823433014,
          'ATC': 0.263743084667417},
    'L': {'TTG': 0.286319859527146,
          'TTA': 0.275534472444779,
          'CTG': 0.110440170850593,
          'CTA': 0.141277445174148,
          'CTT': 0.129115062940288,
          'CTC': 0.0573129890630467},
    'K': {'AAG': 0.423936637198697,
          'AAA': 0.576063362801303},
    'M': {'ATG': 1},
    'F': {'TTT': 0.586126603840976,
          'TTC': 0.413873396159024},
    'P': {'CCG': 0.120626895854398,
          'CCA': 0.417143753704543,
          'CCT': 0.307740315888567,
          'CCC': 0.154489034552491},
    'S': {'AGT': 0.159245398699046,
          'AGC': 0.109749229743856,
          'TCG': 0.0963590866114069,
          'TCA': 0.210157220085731,
          'TCT': 0.264456618519558,
          'TCC': 0.160032446340401},
    'T': {'ACG': 0.135583991997041,
          'ACA': 0.302413913478422,
          'ACT': 0.345237040780705,
          'ACC': 0.216765053743832},
    'W': {'TGG': 1},
    'Y': {'TAT': 0.559573963633711,
          'TAC': 0.440426036366289},
    'V': {'GTG': 0.190897642582249,
          'GTA': 0.208783185960798,
          'GTT': 0.391481704636128,
          'GTC': 0.208837466820824}}

# Complete list of codons.
CODONS = {'AAA': 'K',
          'AAC': 'N',
          'AAG': 'K',
          'AAT': 'N',
          'ACA': 'T',
          'ACC': 'T',
          'ACG': 'T',
          'ACT': 'T',
          'AGA': 'R',
          'AGC': 'S',
          'AGG': 'R',
          'AGT': 'S',
          'ATA': 'I',
          'ATC': 'I',
          'ATG': 'M',
          'ATT': 'I',
          'CAA': 'Q',
          'CAC': 'H',
          'CAG': 'Q',
          'CAT': 'H',
          'CCA': 'P',
          'CCC': 'P',
          'CCG': 'P',
          'CCT': 'P',
          'CGA': 'R',
          'CGC': 'R',
          'CGG': 'R',
          'CGT': 'R',
          'CTA': 'L',
          'CTC': 'L',
          'CTG': 'L',
          'CTT': 'L',
          'GAA': 'E',
          'GAC': 'D',
          'GAG': 'E',
          'GAT': 'D',
          'GCA': 'A',
          'GCC': 'A',
          'GCG': 'A',
          'GCT': 'A',
          'GGA': 'G',
          'GGC': 'G',
          'GGG': 'G',
          'GGT': 'G',
          'GTA': 'V',
          'GTC': 'V',
          'GTG': 'V',
          'GTT': 'V',
          'TAA': '*',
          'TAC': 'Y',
          'TAG': '*',
          'TAT': 'Y',
          'TCA': 'S',
          'TCC': 'S',
          'TCG': 'S',
          'TCT': 'S',
          'TGA': '*',
          'TGC': 'C',
          'TGG': 'W',
          'TGT': 'C',
          'TTA': 'L',
          'TTC': 'F',
          'TTG': 'L',
          'TTT': 'F'}
