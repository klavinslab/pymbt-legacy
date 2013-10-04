'''Read and write DNA sequences.'''
import csv
import os
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import pymbt.sequence
from pymbt.constants import genbank


def NoteNumberError(ValueError):
    pass


def read_dna(path):
    '''Read DNA from file. Uses BioPython and coerces to pymbt format.

    :param path: Full path to input file.
    :type path: str

    '''
    ext = os.path.splitext(path)[1]

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
    dna = pymbt.sequence.DNA(seq.seq.tostring())
    dna.name = seq.name

    # Features
    for feature in seq.features:
        feature_name = feature.qualifiers['label'][0]
        feature_start = int(feature.location.start)
        feature_stop = int(feature.location.end)
        feature_type = _process_feature_type(feature.type)
        if feature.location.strand == -1:
            feature_strand = 1
        else:
            feature_strand = 0
        dna.features.append(pymbt.sequence.Feature(feature_name, feature_start,
                                                   feature_stop, feature_type,
                                                   strand=feature_strand))
    dna.features = sorted(dna.features, key=lambda feature: feature.start)
    try:
        if seq.annotations['data_file_division'] == 'circular':
            dna.topology = 'circular'
        elif seq.annotations['data_file_division'] == 'linear':
            dna.topology = 'linear'
    except KeyError:
        pass

    return dna


def read_sequencing(directory):
    '''Read .seq and .abi/.ab1 results files from a dir.

    :param directory: Path to directory containing sequencing files.
    :type directory: str

    '''
    dirfiles = os.listdir(directory)
    seq_exts = ['.seq', '.abi', '.ab1']
    # Exclude files that aren't sequencing results
    seq_paths = [x for x in dirfiles if os.path.splitext(x)[1] in seq_exts]
    paths = [os.path.join(directory, x) for x in seq_paths]
    sequences = [read_dna(x).set_stranded('ss') for x in paths]

    return sequences


def write_dna(dna, path):
    '''Write DNA to a file (genbank or fasta).

    :param dna: DNA sequence to write to file
    :type dna: pymbt.sequence.DNA
    :param path: file path to write. Has to be genbank or fasta file.
    :type path: str

    '''
    # Check if path filetype is valid, remember for later
    ext = os.path.splitext(path)[1]
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
        bio_strand = 1 if feature.strand == 1 else -1
        location = FeatureLocation(ExactPosition(feature.start),
                                   ExactPosition(feature.stop),
                                   strand=bio_strand)
        ftype = _process_feature_type(feature.feature_type, bio_to_pymbt=False)
        features.append(SeqFeature(location, type=ftype,
                        qualifiers={'label': [feature.name]}))
    # Biopython doesn't like 'None' here
    bio_id = dna.id if dna.id else ''
    # Maximum length of name is 16
    seq = SeqRecord(Seq(str(dna), alphabet=ambiguous_dna), id=bio_id,
                    name=dna.name[0:16], features=features,
                    description=dna.name)

    if filetype == 'genbank':
        SeqIO.write(seq, path, 'genbank')
    elif filetype == 'fasta':
        SeqIO.write(seq, path, 'fasta')


def write_primers(primer_list, path, notes=None):
    """Write a list of primers out to a csv file. The first three columns are
    compatible with the current IDT order form (name, sequence, notes). By
    default there are no notes, which is an optional parameter.

    :param primer_list: A list of primers.
    :type primer_list: pymbt.sequence.Primer list
    :param path: A path to the csv you want to write.
    :type path: str
    :param notes: a list of strings to add for each oligo. Must be same length
                  as primer_list.
    :type notes: str list

    """
    if notes is not None:
        if len(notes) != len(primer_list):
            raise NoteNumberError("Number of primers and notes is not equal.")
    else:
        notes = ["" for x in primer_list]
    with open(path, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["name", "sequence", "notes"])
        for primer, note in zip(primer_list, notes):
            writer.writerow([primer.name, primer.primer(), note])


def _process_feature_type(feature_type, bio_to_pymbt=True):
    '''Translate genbank feature types into usable ones (currently identical).
    The feature table is derived from the official genbank spec (gbrel.txt)
    available at http://www.insdc.org/documents/feature-table

    :param feature_type: feature to convert
    :type feature_type: str
    :param bio_to_pymbt: from pymbt to Biopython (True) or the other direction
                   (False)
    :param bio_to_pymbt: bool

    '''

    err_msg = 'Unrecognized feature type: {}'.format(feature_type)
    if bio_to_pymbt:
        try:
            name = genbank.TO_PYMBT[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    else:
        try:
            name = genbank.TO_BIO[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    return name
