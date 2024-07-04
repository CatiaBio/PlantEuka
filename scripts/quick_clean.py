import gzip
import os
import sys
from collections import Counter

def clean_fasta(fasta_path):
    with gzip.open(fasta_path, 'rt') as file:
        original_content = file.read()
    
    replaced_characters = Counter()
    cleaned_lines = []

    for line in original_content.split('\n'):
        if line.startswith('>'):
            cleaned_lines.append(line)
        else:
            cleaned_line, counts = clean_sequence(line)
            cleaned_lines.append(cleaned_line)
            replaced_characters.update(counts)
    
    cleaned_content = '\n'.join(cleaned_lines)
    
    with gzip.open(fasta_path, 'wt') as file:
        file.write(cleaned_content)

    return replaced_characters

def clean_sequence(sequence):
    replaced_characters = Counter()
    cleaned_sequence = []

    for nucleotide in sequence:
        if nucleotide not in 'ACGTN':
            cleaned_sequence.append('N')
            replaced_characters[nucleotide] += 1
        else:
            cleaned_sequence.append(nucleotide)

    return ''.join(cleaned_sequence), replaced_characters

def process_directory(base_dir):
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.fasta.gz'):
                fasta_path = os.path.join(root, file)
                clean_fasta(fasta_path)
                print(f"Cleaned file: {fasta_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <base_directory>")
        sys.exit(1)

    base_dir = sys.argv[1]
    process_directory(base_dir)
