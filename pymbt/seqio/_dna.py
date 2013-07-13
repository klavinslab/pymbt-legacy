'''Read and write DNA sequences.'''

import os
from Bio import SeqIO
from pymbt import sequence


def read_dna(path):
    '''Read DNA from file. Uses BioPython and coerces to pymbt format.

    :param path: Full path to input file.
    :type path: str

    '''
    base, ext = os.path.splitext(path)

    genbank_exts = ['.gb', '.ape']
    fasta_exts = ['.fasta', '.fa', '.seq']
    abi_exts = ['.abi', '.ab1']

    if any(ext == extension for extension in genbank_exts):
        file_format = 'genbank'
    elif any(ext == extension for extension in fasta_exts):
        file_format = 'fasta'
    elif any(ext == extension for extension in abi_exts):
        file_format = 'abi'
    else:
        raise ValueError('File format not recognized.')

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
    '''Read .seq and .abi/.ab1 results files from a dir.

    :param dirpath: Path to directory containing sequencing files.
    :type dirpath: str

    '''
    seq_exts = ['.seq', '.abi', '.ab1']
    dirfiles = os.listdir(dirpath)
    seq_paths = [x for x in dirfiles if os.path.splitext(x)[1] in seq_exts]
    sequences = [read_dna(dirpath + x).set_stranded('ss') for x in seq_paths]

    return sequences


def _process_feature_type(feature_type):
    '''Translate BioPython / genbank feature types into usable ones.

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
