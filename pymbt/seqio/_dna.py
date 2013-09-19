'''Read and write DNA sequences.'''
import os
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import pymbt.sequence


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
        feature_strand = max(feature.location.strand, 0)
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
    to_pymbt = {'-': '-',
                '-10_signal': '-10_signal',
                '-35_signal': '-35_signal',
                "3'UTR": "3'UTR",
                "3'clip": "3'clip",
                "5'UTR": "5'UTR",
                "5'clip": "5'clip",
                'CAAT_signal': 'CAAT_signal',
                'CDS': 'CDS',
                'C_region': 'C_region',
                'D-loop': 'D-loop',
                'D_segment': 'D_segment',
                'GC_signal': 'GC_signal',
                'J_region': 'J_region',
                'LTR': 'LTR',
                'N_region': 'N_region',
                'RBS': 'RBS',
                'STS': 'STS',
                'S_region': 'S_region',
                'TATA_signal': 'TATA_signal',
                'V_region': 'V_region',
                'allele': 'allele',
                'attenuator': 'attenuator',
                'conflict': 'conflict',
                'enhancer': 'enhancer',
                'exon': 'exon',
                'gene': 'gene',
                'iDNA': 'iDNA',
                'intron': 'intron',
                'mRNA': 'mRNA',
                'mat_peptide': 'mat_peptide',
                'misc_RNA': 'misc_RNA',
                'misc_binding': 'misc_binding',
                'misc_difference': 'misc_difference',
                'misc_feature': 'misc_feature',
                'misc_recomb': 'misc_recomb',
                'misc_signal': 'misc_signal',
                'misc_structure': 'misc_structure',
                'modified_base': 'modified_base',
                'mutation ': 'mutation ',
                'old_sequence': 'old_sequence',
                'polyA_signal': 'polyA_signal',
                'polyA_site': 'polyA_site',
                'precursor_RNA': 'precursor_RNA',
                'prim_transcript': 'prim_transcript',
                'primer': 'primer',
                'primer_bind': 'primer_bind',
                'promoter': 'promoter',
                'protein_bind': 'protein_bind',
                'rRNA': 'rRNA',
                'rep_origin': 'rep_origin',
                'repeat_region': 'repeat_region',
                'repeat_unit': 'repeat_unit',
                'satellite': 'satellite',
                'scRNA': 'scRNA',
                'sig_peptide': 'sig_peptide',
                'snRNA': 'snRNA',
                'source': 'source',
                'stem_loop': 'stem_loop',
                'tRNA ': 'tRNA ',
                'terminator': 'terminator',
                'transit_peptide': 'transit_peptide',
                'transposon': 'transposon',
                'unsure': 'unsure',
                'variation ': 'variation '}
    to_bio = {value: key for key, value in to_pymbt}

    err_msg = 'Unrecognized feature type: {}'.format(feature_type)
    if bio_to_pymbt:
        try:
            name = to_pymbt[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    else:
        try:
            name = to_bio[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    return name
