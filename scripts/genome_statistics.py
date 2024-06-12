#!/usr/bin/env python3

"""
Description:
This script calculates nucleotide statistics for each .fasta.gz file in a given directory and outputs the results to TSV files.
It processes each file to determine counts of nucleotides A, T, C, G, N, total length, and percentage content of A:T and C:G pairs.

Example:
python3 scripts/genome_statistics.py genomes/2405_chloroplast_merged results/2405_chloroplast
"""

# Libraries 
import sys
import gzip
import os
from collections import Counter

def ensure_directory_exists(directory):
    """Ensure the specified directory exists, and if not, create it."""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")

def parse_fasta(fasta_file):
    """
    Parse a FASTA file compressed with gzip.

    Args:
    fasta_file (str): Path to the gzipped FASTA file.

    Yields:
    tuple: Tuple containing the header and sequence of each entry in the FASTA file.
    """
    header = None
    sequence = []
    with gzip.open(fasta_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, ''.join(sequence)
                header = line.split()[0][1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)

def calculate_nucleotide_statistics(fasta_file, output_file, taxonomy, name):
    """
    Calculate nucleotide statistics for sequences in a FASTA file and write the results to a TSV file.

    Args:
    fasta_file (str): Path to the gzipped FASTA file to process.
    output_file (str): Path to the TSV output file where statistics will be saved.
    taxonomy (str): Taxonomy extracted from the filename.
    name (str): Name extracted from the filename.
    """
    total_counts = Counter()
    num_sequences = 0
    with open(output_file, 'w') as out:
        out.write("taxonomy\tname\taccession.number\tA\tT\tC\tG\tN\tlength\tA:T%\tC:G%\n")
        for header, sequence in parse_fasta(fasta_file):
            counts = Counter(sequence.upper())
            length = len(sequence)
            a, t, c, g, n = counts['A'], counts['T'], counts['C'], counts['G'], counts['N']
            at_content = ((a + t) / length * 100) if length else 0
            cg_content = ((c + g) / length * 100) if length else 0
            out.write(f"{taxonomy}\t{name}\t{header}\t{a}\t{t}\t{c}\t{g}\t{n}\t{length}\t{at_content:.2f}\t{cg_content:.2f}\n")
            total_counts += counts
            num_sequences += 1
        
def process_directory(directory, output):
    """Process all .fasta.gz files in the directory to calculate statistics."""
    ensure_directory_exists(output)  # Ensure output directory exists
    for file in os.listdir(directory):
        if file.endswith(".fasta.gz"):
            fasta_file = os.path.join(directory, file)
            output_file = os.path.join(output, f"{file.rsplit('.fasta.gz', 1)[0]}_stats.tsv")
            # Extract taxonomy and name from the filename
            taxonomy, name = file.rsplit('.', 2)[0].split('_', 1)
            calculate_nucleotide_statistics(fasta_file, output_file, taxonomy, name)
            print(f"Processed: {fasta_file} -> {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fasta_folder_stats.py <input_folder> <output_folder>")
        sys.exit(1)
    organelle = sys.argv[1]    
    input_folder = f"{organelle}/genomes/merged"
    output_folder = f"{organelle}/results/stats"
    process_directory(input_folder, output_folder)