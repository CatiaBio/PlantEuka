#!/usr/bin/env python3

import os
import gzip
from collections import Counter
import csv
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 2:
    print("Usage: python3 scripts/sort_genomes.py <organelle>")
    sys.exit(1)

organelle = sys.argv[1]

# Paths
genomes_folder = f"{organelle}/genomes/original"
accession_taxid_file = f"{organelle}/other/accessions_taxid.txt"
lineage_file = "other/lineage.tsv"
output_file = f"{organelle}/other/data_full_info.tsv"

# Read accession_taxid data
accession_taxid = {}
with open(accession_taxid_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            accession_taxid[parts[0]] = parts[1]
        else:
            print(f"Skipping line in accession_taxid_file: {line.strip()}")

# Read lineage data
lineage_data = {}
with open(lineage_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    headers = next(reader)
    for row in reader:
        if len(row) > 1:
            lineage_data[row[0]] = row[1:]
        else:
            print(f"Skipping line in lineage_file: {row}")

# Function to count nucleotides in a sequence
def count_nucleotides(sequence):
    counts = Counter(sequence)
    total_length = sum(counts.values())
    a_count = counts.get('A', 0)
    t_count = counts.get('T', 0)
    c_count = counts.get('C', 0)
    g_count = counts.get('G', 0)
    n_count = total_length - (a_count + t_count + c_count + g_count)  # Count non-ACGT as N
    at_content = ((a_count + t_count) / total_length * 100) if total_length else 0
    cg_content = ((c_count + g_count) / total_length * 100) if total_length else 0
    return a_count, t_count, c_count, g_count, n_count, total_length, at_content, cg_content

# Process fasta files and gather information
output_data = []
for filename in os.listdir(genomes_folder):
    if filename.endswith(".fasta.gz"):
        accession_version = filename.split('.')[0]
        filepath = os.path.join(genomes_folder, filename)
        
        try:
            # Read fasta file
            with gzip.open(filepath, 'rt') as f:
                lines = f.readlines()
                if lines:
                    header = lines[0].strip()
                    accession_in_header = header.split(' ')[0].replace('>', '')
                    species = ' '.join(header.split(' ')[1:3])  # Assuming species name is the second and third word
                    sequence = ''.join([line.strip() for line in lines[1:]])
                    
                    # Count nucleotides
                    a_count, t_count, c_count, g_count, n_count, total_length, at_content, cg_content = count_nucleotides(sequence)
                    
                    # Get taxid
                    taxid = accession_taxid.get(accession_in_header, 'NA')
                    
                    # Get lineage information
                    lineage_info = lineage_data.get(taxid, ['NA'] * 6)
                    
                    # Append to output data
                    output_data.append([
                        accession_version, taxid, *lineage_info, 
                        a_count, t_count, c_count, g_count, n_count, total_length, at_content, cg_content
                    ])
                else:
                    print(f"No content in file: {filepath}")
        except Exception as e:
            print(f"Error processing file {filepath}: {e}")

# Write output to TSV file
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    # Write header
    writer.writerow([
        'accession.version', 'taxid', 'species', 'genus', 'family', 'order', 'class', 'phylum', 
        'A', 'T', 'C', 'G', 'N', 'length', 'A:T', 'C:G'
    ])
    # Write data
    writer.writerows(output_data)

print(f"Output written to {output_file}")
