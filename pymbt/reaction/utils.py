'''
Utilities for reactions.

'''

from pymbt import data
from pymbt import sequence


def convert_sequence(seq, from_material, to_material):
    '''
    Translate a DNA sequence into peptide sequence.

    :param seq: DNA or RNA sequence.
    :type seq: DNA or RNA
    :param from_material: material to convert ('rna', or 'dna')
    :type from_material: str
    :param to_material: material to which to convert ('rna', 'dna', or
                        'peptide').
    :type to_material: str

    '''

    if from_material == 'dna' and to_material == 'rna':
        # Can't transcribe a gap
        if sequence.DNA('-') in seq:
            raise ValueError('Cannot transcribe gapped DNA')
        # Convert DNA chars to RNA chars
        origin = data.common.ALPHABETS['dna'][:-1]
        destination = data.common.ALPHABETS['rna']
        code = dict(zip(origin, destination))
        converted = ''.join(code.get(str(k), str(k)) for k in seq)
        # Instantiate RNA object
        converted = sequence.RNA(converted)
    elif from_material == 'rna' and to_material == 'dna':
        # Convert RNA chars to DNA chars
        origin = data.common.ALPHABETS['rna']
        destination = data.common.ALPHABETS['dna'][:-1]
        code = dict(zip(origin, destination))
        converted = ''.join(code.get(str(k), str(k)) for k in seq)
        # Instantiate DNA object
        converted = sequence.DNA(converted)
    elif from_material == 'rna' and to_material == 'peptide':
        # Make a list for easier processing
        seq_list = list(str(seq))

        # Convert to peptide until stop codon is found.
        converted = []
        while True:
            if len(seq_list) >= 3:
                base_1 = seq_list.pop(0)
                base_2 = seq_list.pop(0)
                base_3 = seq_list.pop(0)
                codon = ''.join(base_1 + base_2 + base_3).upper()
                amino_acid = data.common.CODONS[codon]
                # Stop when stop codon is found
                if amino_acid == '*':
                    break
                converted.append(amino_acid)
            else:
                break
        converted = ''.join(converted)
        converted = sequence.Peptide(converted)
    else:
        msg1 = 'Conversion from '
        msg2 = '{0} to {1} is not supported.'.format(from_material,
                                                     to_material)
        raise ValueError(msg1 + msg2)

    return converted
