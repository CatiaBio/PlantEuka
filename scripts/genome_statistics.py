#!/usr/bin/env python3

import sys
import gzip
import os
from collections import Counter

def parse_fasta(fasta_file):
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

def calculate_nucleotide_statistics(fasta_file, output_file):
    total_counts = Counter()
    num_sequences = 0
    with open(output_file, 'w') as out:
        out.write("id\tA\tT\tC\tG\tN\tlength\tA:T %\tC:G %\n")
        for header, sequence in parse_fasta(fasta_file):
            counts = Counter(sequence.upper())
            length = len(sequence)
            a, t, c, g, n = counts['A'], counts['T'], counts['C'], counts['G'], counts['N']
            at_content = ((a + t) / length * 100) if length else 0
            cg_content = ((c + g) / length * 100) if length else 0
            out.write(f"{header}\t{a}\t{t}\t{c}\t{g}\t{n}\t{length}\t{at_content:.2f}\t{cg_content:.2f}\n")
            total_counts += counts
            num_sequences += 1
        
        if num_sequences > 0:
            avg_a = total_counts['A'] / num_sequences
            avg_t = total_counts['T'] / num_sequences
            avg_c = total_counts['C'] / num_sequences
            avg_g = total_counts['G'] / num_sequences
            avg_n = total_counts['N'] / num_sequences
            avg_length = sum(total_counts.values()) / num_sequences
            avg_at_content = ((avg_a + avg_t) / avg_length * 100) if avg_length else 0
            avg_cg_content = ((avg_c + avg_g) / avg_length * 100) if avg_length else 0
            out.write(f"Average\t{avg_a:.2f}\t{avg_t:.2f}\t{avg_c:.2f}\t{avg_g:.2f}\t{avg_n:.2f}\t{avg_length:.2f}\t{avg_at_content:.2f}\t{avg_cg_content:.2f}\n")

def process_directory(directory):
    for file in os.listdir(directory):
        if file.endswith(".fasta.gz"):
            fasta_file = os.path.join(directory, file)
            output_file = os.path.join(directory, f"{file.rsplit('.fa', 1)[0]}_stats.tsv")
            calculate_nucleotide_statistics(fasta_file, output_file)
            print(f"Processed: {fasta_file} -> {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fasta_folder_stats.py <input_folder>")
        sys.exit(1)
    input_folder = sys.argv[1]
    process_directory(input_folder)
