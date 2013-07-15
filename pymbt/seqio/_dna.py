'''Read and write DNA sequences.'''

import os
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
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


def write_dna(dna, path):
    '''Write DNA to a file (genbank or fasta).

    :param dna: DNA sequence to write to file
    :type dna: pymbt.sequence.DNA
    :param path: file path to write. Has to be genbank or fasta file.
    :type path: str

    '''
    # Check if path filetype is valid, remember for later
    base, ext = os.path.splitext(path)
    if ext == '.gb' or ext == '.ape':
        filetype = 'genbank'
    elif ext == '.fa' or ext == '.fasta':
        filetype = 'fasta'
    else:
        raise ValueError('Only genbank or fasta files are supported.')

    # Convert features to Biopython form
    # Information lost on conversion:
    #     specificity of feature typ
    #     strandedness
    #     topology
    features = []
    for feature in dna.features:
        location = FeatureLocation(ExactPosition(feature.start),
                                   ExactPosition(feature.stop),
                                   strand=feature.strand)
        ftype = _process_feature_type(feature.feature_type, bio_to_pymbt=False)
        features.append(SeqFeature(location, type=ftype,
                        qualifiers={'label': [feature.name]}))
    # Biopython doesn't like 'None' here
    bio_id = dna.id if dna.id else ''
    # Maximum length of name is 16
    seq = SeqRecord(Seq(dna.top, alphabet=ambiguous_dna), id=bio_id,
                    name=dna.name[0:16], features=features,
                    description=dna.name)

    if filetype == 'genbank':
        SeqIO.write(seq, path, 'genbank')
    elif filetype == 'fasta':
        SeqIO.write(seq, path, 'fasta')


def _process_feature_type(feature_type, bio_to_pymbt=True):
    '''Translate BioPython / genbank feature types into usable ones.

    :param feature_type: feature to convert
    :type feature_type: str
    :param bio_to_pymbt: from pymbt to Biopython (True) or the other direction
                   (False)
    :param bio_to_pymbt: bool

    '''
    to_bio = {'misc_feature': 'misc',
              'CDS': 'coding',
              'gene': 'coding',
              'site': 'misc',
              'primer_bind': 'primer',
              'rep_origin': 'origin',
              'promoter': 'promoter',
              'terminator': 'terminator'}
    to_pymbt = {'misc': 'misc_feature',
                'coding': 'CDS',
                'primer': 'primer_bind',
                'origin': 'rep_origin',
                'promoter': 'promoter',
                'terminator': 'terminator'}

    err_msg = 'Unrecognized feature type: {}'.format(feature_type)
    if bio_to_pymbt:
        try:
            name = to_bio[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    else:
        try:
            name = to_pymbt[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    return name
