'''Defines data and parameters in an easily resuable format.'''

# Common sequence alphabets.
ALPHABETS = {
    'dna': 'ATGCNatgcn-',
    'rna': 'AUGCNaugcn',
    'peptide': 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy'}

COMPLEMENTS = {
    'dna': ('ATGCNatgcn-', 'TACGNtacgn-'),
    'rna': ('AUGCNaugcn', 'UACGNuacgn')}


# The standard codon table.
CODON_TABLE = {
    'A': ['GCG', 'GCA', 'GCU', 'GCC'],
    'R': ['AGG', 'AGA', 'CGG', 'CGA', 'CGU', 'CGC'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    '*': ['UGA', 'UAG', 'UAA'],
    'Q': ['CAG', 'CAA'],
    'E': ['GAG', 'GAA'],
    'G': ['GGG', 'GGA', 'GGU', 'GGC'],
    'H': ['CAU', 'CAC'],
    'I': ['AUA', 'AUU', 'AUC'],
    'L': ['UUG', 'UUA', 'CUG', 'CUA', 'CUU', 'CUC'],
    'K': ['AAG', 'AAA'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCG', 'CCA', 'CCU', 'CCC'],
    'S': ['AGU', 'AGC', 'UCG', 'UCA', 'UCU', 'UCC'],
    'T': ['ACG', 'ACA', 'ACU', 'ACC'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUG', 'GUA', 'GUU', 'GUC']}


# Saccharomyces cerevisiae
# source: http://www.kazusa.or.jp/codon/
# (which cites GenBank, i.e. yeast genome project CDS database)
CODON_FREQ = {
    'sc': {
        'GCG': 0.109972396541529,
        'GCA': 0.288596474496094,
        'GCU': 0.377014739102356,
        'GCC': 0.224416389860021,
        'AGG': 0.208564104515562,
        'AGA': 0.481137590939125,
        'CGG': 0.0392677130215486,
        'CGA': 0.0676728924436203,
        'CGU': 0.144572019635586,
        'CGC': 0.0587856794445578,
        'AAU': 0.589705127199784,
        'AAC': 0.410294872800217,
        'GAU': 0.65037901553924,
        'GAC': 0.34962098446076,
        'UGU': 0.629812614586062,
        'UGC': 0.370187385413938,
        'UGA': 0.303094329334787,
        'UAG': 0.225736095965104,
        'UAA': 0.471169574700109,
        'CAG': 0.307418833439535,
        'CAA': 0.692581166560465,
        'GAG': 0.296739610207218,
        'GAA': 0.703260389792782,
        'GGG': 0.119057918187951,
        'GGA': 0.215422869017838,
        'GGU': 0.472217600813099,
        'GGC': 0.193301611981112,
        'CAU': 0.636710255236351,
        'CAC': 0.363289744763649,
        'AUA': 0.273331091899568,
        'AUU': 0.462925823433014,
        'AUC': 0.263743084667417,
        'UUG': 0.286319859527146,
        'UUA': 0.275534472444779,
        'CUG': 0.110440170850593,
        'CUA': 0.141277445174148,
        'CUU': 0.129115062940288,
        'CUC': 0.0573129890630467,
        'AAG': 0.423936637198697,
        'AAA': 0.576063362801303,
        'AUG': 1,
        'UUU': 0.586126603840976,
        'UUC': 0.413873396159024,
        'CCG': 0.120626895854398,
        'CCA': 0.417143753704543,
        'CCU': 0.307740315888567,
        'CCC': 0.154489034552491,
        'AGU': 0.159245398699046,
        'AGC': 0.109749229743856,
        'UCG': 0.0963590866114069,
        'UCA': 0.210157220085731,
        'UCU': 0.264456618519558,
        'UCC': 0.160032446340401,
        'ACG': 0.135583991997041,
        'ACA': 0.302413913478422,
        'ACU': 0.345237040780705,
        'ACC': 0.216765053743832,
        'UGG': 1,
        'UAU': 0.559573963633711,
        'UAC': 0.440426036366289,
        'GUG': 0.190897642582249,
        'GUA': 0.208783185960798,
        'GUU': 0.391481704636128,
        'GUC': 0.208837466820824}}

# Codon usage for S cerevisiae
CODON_FREQ_SC_NESTED = {
    'A': {'GCG': 0.109972396541529,
          'GCA': 0.288596474496094,
          'GCU': 0.377014739102356,
          'GCC': 0.224416389860021},
    'R': {'AGG': 0.208564104515562,
          'AGA': 0.481137590939125,
          'CGG': 0.0392677130215486,
          'CGA': 0.0676728924436203,
          'CGU': 0.144572019635586,
          'CGC': 0.0587856794445578},
    'N': {'AAU': 0.589705127199784,
          'AAC': 0.410294872800217},
    'D': {'GAU': 0.65037901553924,
          'GAC': 0.34962098446076},
    'C': {'UGU': 0.629812614586062,
          'UGC': 0.370187385413938},
    '*': {'UGA': 0.303094329334787,
          'UAG': 0.225736095965104,
          'UAA': 0.471169574700109},
    'Q': {'CAG': 0.307418833439535,
          'CAA': 0.692581166560465},
    'E': {'GAG': 0.296739610207218,
          'GAA': 0.703260389792782},
    'G': {'GGG': 0.119057918187951,
          'GGA': 0.215422869017838,
          'GGU': 0.472217600813099,
          'GGC': 0.193301611981112},
    'H': {'CAU': 0.636710255236351,
          'CAC': 0.363289744763649},
    'I': {'AUA': 0.273331091899568,
          'AUU': 0.462925823433014,
          'AUC': 0.263743084667417},
    'L': {'UUG': 0.286319859527146,
          'UUA': 0.275534472444779,
          'CUG': 0.110440170850593,
          'CUA': 0.141277445174148,
          'CUU': 0.129115062940288,
          'CUC': 0.0573129890630467},
    'K': {'AAG': 0.423936637198697,
          'AAA': 0.576063362801303},
    'M': {'AUG': 1},
    'F': {'UUU': 0.586126603840976,
          'UUC': 0.413873396159024},
    'P': {'CCG': 0.120626895854398,
          'CCA': 0.417143753704543,
          'CCU': 0.307740315888567,
          'CCC': 0.154489034552491},
    'S': {'AGU': 0.159245398699046,
          'AGC': 0.109749229743856,
          'UCG': 0.0963590866114069,
          'UCA': 0.210157220085731,
          'UCU': 0.264456618519558,
          'UCC': 0.160032446340401},
    'T': {'ACG': 0.135583991997041,
          'ACA': 0.302413913478422,
          'ACU': 0.345237040780705,
          'ACC': 0.216765053743832},
    'W': {'UGG': 1},
    'Y': {'UAU': 0.559573963633711,
          'UAC': 0.440426036366289},
    'V': {'GUG': 0.190897642582249,
          'GUA': 0.208783185960798,
          'GUU': 0.391481704636128,
          'GUC': 0.208837466820824}}

# Complete list of codons.
CODONS = {'AAA': 'K',
          'AAC': 'N',
          'AAG': 'K',
          'AAU': 'N',
          'ACA': 'T',
          'ACC': 'T',
          'ACG': 'T',
          'ACU': 'T',
          'AGA': 'R',
          'AGC': 'S',
          'AGG': 'R',
          'AGU': 'S',
          'AUA': 'I',
          'AUC': 'I',
          'AUG': 'M',
          'AUU': 'I',
          'CAA': 'Q',
          'CAC': 'H',
          'CAG': 'Q',
          'CAU': 'H',
          'CCA': 'P',
          'CCC': 'P',
          'CCG': 'P',
          'CCU': 'P',
          'CGA': 'R',
          'CGC': 'R',
          'CGG': 'R',
          'CGU': 'R',
          'CUA': 'L',
          'CUC': 'L',
          'CUG': 'L',
          'CUU': 'L',
          'GAA': 'E',
          'GAC': 'D',
          'GAG': 'E',
          'GAU': 'D',
          'GCA': 'A',
          'GCC': 'A',
          'GCG': 'A',
          'GCU': 'A',
          'GGA': 'G',
          'GGC': 'G',
          'GGG': 'G',
          'GGU': 'G',
          'GUA': 'V',
          'GUC': 'V',
          'GUG': 'V',
          'GUU': 'V',
          'UAA': '*',
          'UAC': 'Y',
          'UAG': '*',
          'UAU': 'Y',
          'UCA': 'S',
          'UCC': 'S',
          'UCG': 'S',
          'UCU': 'S',
          'UGA': '*',
          'UGC': 'C',
          'UGG': 'W',
          'UGU': 'C',
          'UUA': 'L',
          'UUC': 'F',
          'UUG': 'L',
          'UUU': 'F'}
