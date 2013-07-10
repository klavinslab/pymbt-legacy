'''
Read and write DNA sequences.

'''

from os import listdir
from Bio import SeqIO
from pymbt import sequence


def read_dna(path, file_format):
    '''
    Read DNA from file. Uses BioPython's tools and coerces to pymbt format.

    :param path: Full path to input file.
    :type path: str
    :param file_format: BioPython-compatible format string, e.g. 'fasta',
                        'genbank'.
    :type file_format: str
    '''

    seq = SeqIO.read(path, file_format)
    dna = sequence.DNA(seq.seq.tostring())
    dna.name = seq.name

    # Features
    for feature in seq.features:
        feature_name = feature.qualifiers['label'][0]
        feature_start = int(feature.location.start)
        feature_stop = int(feature.location.end)
        feature_type = _process_feature_type(feature.type)
        feature_strand = max(feature.location.strand, 0)
        dna.features.append(sequence.Feature(feature_name, feature_start,
                                             feature_stop, feature_type,
                                             strand=feature_strand))
    dna.features = sorted(dna.features, key=lambda feature: feature.start)

    return dna


def read_sequencing(dirpath):
    '''
    Read .seq and .abi/.ab1 results files from a dir.

    :param dirpath: Path to directory containing sequencing files.
    :type dirpath: str

    '''

    seq_paths = [x for x in listdir(dirpath) if x.endswith('.seq')]
    abi_paths = [x for x in listdir(dirpath) if x.endswith('.abi')]
    abi_paths += [x for x in listdir(dirpath) if x.endswith('.ab1')]
    seq_seqs = [read_dna(dirpath + x, 'fasta') for x in seq_paths]
    abi_seqs = [read_dna(dirpath + x, 'abi') for x in abi_paths]
    sequences = seq_seqs + abi_seqs
    sequences = [seq.set_stranded('ss') for seq in sequences]

    return sequences


def _process_feature_type(feature_type):
    '''
    Translate BioPython / genbank feature types into those used by
    pymbt.sequence.Feature

    '''

    conversion = {'misc_feature': 'misc',
                  'CDS': 'coding',
                  'gene': 'coding',
                  'site': 'misc',
                  'primer_bind': 'primer',
                  'rep_origin': 'origin',
                  'promoter': 'promoter',
                  'terminator': 'terminator'}
    try:
        name = conversion[feature_type]
    except KeyError:
        raise ValueError('Unrecognized feature type: {}'.format(feature_type))
    return name
