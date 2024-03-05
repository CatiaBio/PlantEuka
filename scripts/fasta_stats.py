#!/usr/bin/env python3

import sys
import gzip
from collections import Counter

def parse_fasta(fasta_file):
    """Parse a gzipped FASTA file and yield sequence headers and sequences."""
    header = None
    sequence = []
    # Determine if the file is gzipped based on the extension and choose the appropriate open function
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    with open_func(fasta_file, 'rt') as f:  # 'rt' mode for reading text from gzip
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, ''.join(sequence)
                # Extract the specific part of the header (e.g., "NC_023337.1")
                header_part = line.split()[0][1:]  # Split by space and remove '>'
                header = header_part
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)

def calculate_nucleotide_statistics(fasta_file, output_file):
    """Calculate nucleotide statistics for each sequence in a gzipped FASTA file."""
    with open(output_file, 'w') as out:
        # Write the header row
        out.write("id\tA\tT\tC\tG\tN\tlength\n")
        for header, sequence in parse_fasta(fasta_file):
            # Count the occurrences of each nucleotide
            counts = Counter(sequence.upper())
            # Calculate the length of the sequence
            length = len(sequence)
            # Write the statistics to the output file
            out.write(f"{header}\t{counts['A']}\t{counts['T']}\t{counts['C']}\t{counts['G']}\t{counts['N']}\t{length}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fasta_stats.py <input_fasta_file> <output_tsv_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    calculate_nucleotide_statistics(fasta_file, output_file)

