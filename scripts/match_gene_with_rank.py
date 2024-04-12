#!/usr/bin/env python3

"""
Description:
This script matches gene information from specified files with rank data and saves the results.
The gene's organelle type can be specified, allowing for dynamic handling of different organelles.

Usage:
./match_gene_with_rank.py <gene_name> <organelle_type>

Arguments:
<gene_name>: Name of the gene for which data is processed (e.g., rbcL).
<organelle_type>: Type of organelle (e.g., chloroplast, mitochondrion).

Example usage following PlantEuka folder organization:
./scripts/match_gene_with_rank.py rbcL chloroplast
"""

# Libraries 
import csv
import glob
import os
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) < 3:
    print("Usage: ./match_gene_with_rank.py <gene_name> <organelle_type>")
    sys.exit(1)  

# Extract command-line arguments
gene_name = sys.argv[1]
organelle_type = sys.argv[2]

# Path settings
base_dir = f"genomes/{organelle_type}/merged"
test_gene_file_path = os.path.join(base_dir, f"{gene_name}/{gene_name}_gene_info.tsv")
input_dir_path = os.path.join(base_dir, 'acc_numb')
results_base_dir = os.path.join(base_dir, gene_name)
missing_files_dir = os.path.join(results_base_dir, 'missing')

# Ensure the directories exist
os.makedirs(results_base_dir, exist_ok=True)
os.makedirs(missing_files_dir, exist_ok=True)

# List all .tsv files in the directory
tsv_files = glob.glob(os.path.join(input_dir_path, '*.tsv'))

def process_file(accession_file_path, test_gene_file_path):
    """
    Process each TSV file to match gene information with accession numbers and save results.

    Args:
    accession_file_path (str): Path to the TSV file containing accession numbers.
    test_gene_file_path (str): Path to the gene info file to match against the accession numbers.
    """
    base_name = os.path.basename(accession_file_path).replace('.tsv', '')
    output_file_path = os.path.join(results_base_dir, f"{base_name}_{gene_name}_gene.tsv")
    missing_file_path = os.path.join(missing_files_dir, f"{base_name}_missing.tsv")
    
    accession_numbers = set()
    with open(accession_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            accession_numbers.add(row[0].split(' ')[0])

    matching_entries = []
    found_accessions = set()

    with open(test_gene_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            acc_version, location_info = row[0].split(' ', 1)
            strand = location_info.split(']')[-1]
            if acc_version in accession_numbers:
                start_end = location_info.split(']')[0].strip('[')
                start, end = start_end.split(':')
                matching_entries.append({'acc_version': acc_version, 'start_pos': int(start) + 1, 'end_pos': int(end), 'orientation': strand})
                found_accessions.add(acc_version)

    missing_accessions = accession_numbers - found_accessions

    if not missing_accessions:
        with open(output_file_path, 'w', newline='') as output_file:
            fieldnames = ['acc_version', 'start_pos', 'end_pos', 'orientation']
            writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for entry in matching_entries:
                writer.writerow(entry)
    else:
        with open(missing_file_path, 'w', newline='') as missing_file:
            writer = csv.writer(missing_file, delimiter='\t')
            for acc in missing_accessions:
                writer.writerow([acc])
        print(f"Created missing accession numbers file: {missing_file_path}")

# Process each .tsv file in the directory
for tsv_file in tsv_files:
    if f'_{gene_name}_gene.tsv' not in tsv_file:
        process_file(tsv_file, test_gene_file_path)