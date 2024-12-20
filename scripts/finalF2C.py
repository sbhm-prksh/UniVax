import sys
import pandas as pd
import re
import os

def is_natural_sequence(sequence):
    """Strictly check if sequence contains only standard amino acids"""
    standard_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in standard_aa for aa in sequence)

def parse_fasta(file):
    data = {
        'Protein ID': [],
        'Allergen Test': [],
        'Antigen Test': [],
        'Signal P': [],
        'Sequence': []
    }
    total_sequences = 0
    natural_sequences = 0

    with open(file, 'r') as f:
        protein_id = ""
        sequence = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if protein_id:
                    total_sequences += 1
                    full_sequence = ''.join(sequence).upper()
                    if is_natural_sequence(full_sequence):
                        natural_sequences += 1
                        data['Protein ID'].append(protein_id)
                        data['Sequence'].append(full_sequence)
                        data['Allergen Test'].append("Pending")
                        data['Antigen Test'].append("Pending")
                        data['Signal P'].append("Pending")
                protein_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

        if protein_id:
            total_sequences += 1
            full_sequence = ''.join(sequence).upper()
            if is_natural_sequence(full_sequence):
                natural_sequences += 1
                data['Protein ID'].append(protein_id)
                data['Sequence'].append(full_sequence)
                data['Allergen Test'].append("Pending")
                data['Antigen Test'].append("Pending")
                data['Signal P'].append("Pending")

    return data, total_sequences, natural_sequences

def fasta_to_csv(fasta_file, csv_file):
    data, total_seqs, natural_seqs = parse_fasta(fasta_file)
    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False)

    print(f"Total sequences: {total_seqs}")
    print(f"Natural sequences: {natural_seqs}")
    print(f"Discarded sequences: {total_seqs - natural_seqs}")
    print(f"Success: Converted {fasta_file} -> {csv_file}")

if __name__ == "__main__":
    # Read arguments from the command line
    if len(sys.argv) != 3:
        print("Usage: python finalF2C.py <input_fasta> <output_csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist!")
        sys.exit(1)

    fasta_to_csv(input_file, output_file)
