#!/usr/bin/env python3

"""
Description:
This script modifies gene sequences based on orientation and start position from gene info.
The script dynamically handles different organelles and gene names as specified by the user.

Usage:
python script_name.py <gene_name> <organelle_type>

Arguments:
<gene_name>: The name of the gene to process (e.g., rbcL).
<organelle_type>: The organelle type, e.g., chloroplast or mitochondrion.

Example:
python script_name.py rbcL chloroplast
"""

# Libraries 
from Bio.Seq import Seq
from Bio import SeqIO
import csv
import gzip
import glob
import os
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) < 3:
    print("Usage: python script_name.py <gene_name> <organelle_type>")
    sys.exit(1)

# Extract command-line arguments
gene_name = sys.argv[1]  
organelle_type = sys.argv[2]  

# Dynamically set paths based on organelle type
fasta_dir_path = f'genomes/{organelle_type}/merged'
tsv_dir_path = f'genomes/{organelle_type}/merged/{gene_name}'
modified_dir_path = f'genomes/{organelle_type}/merged/modified_{gene_name}'

os.makedirs(modified_dir_path, exist_ok=True)  # Ensure the modified directory exists

def load_gene_info(gene_info_path):
    """
    Load gene information from a TSV file.

    Args:
    gene_info_path (str): The file path to the gene information TSV file.

    Returns:
    dict: A dictionary with accession versions as keys and tuples of start position, end position, and orientation as values.
    """
    gene_info = {}
    with open(gene_info_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            acc_version = row['acc_version']
            start_pos = int(row['start_pos'])
            end_pos = int(row['end_pos'])
            orientation = row['orientation']
            gene_info[acc_version] = (start_pos, end_pos, orientation)
    return gene_info

def modify_sequences(fasta_file_path, gene_info, output_fasta_path):
    """
    Modify sequences by applying reverse complement and rotating to start at the gene's position.

    Args:
    fasta_file_path (str): The path to the input FASTA file.
    gene_info (dict): A dictionary containing gene start and end positions and orientation.
    output_fasta_path (str): The path to the output modified FASTA file.
    """
    with gzip.open(fasta_file_path, 'rt') as file, gzip.open(output_fasta_path, 'wt') as output_file:
        for record in SeqIO.parse(file, 'fasta'):
            acc_version = record.id.split(' ')[0]  # Use the full accession number including version
            if acc_version in gene_info:
                start_pos, end_pos, orientation = gene_info[acc_version]

                # Apply reverse complement first if orientation is negative
                if orientation == '(-)':
                    record.seq = record.seq.reverse_complement()
                    # Adjust start and end positions for the new sequence orientation
                    start_pos, end_pos = len(record.seq) - end_pos + 1, len(record.seq) - start_pos + 1

                # Now, rotate the sequence to start at the gene's start position
                rotated_seq = record.seq[start_pos-1:] + record.seq[:start_pos-1]
                record.seq = rotated_seq

                SeqIO.write([record], output_file, 'fasta')
            else:
                print(f"Gene info not found for {acc_version} in file {os.path.basename(fasta_file_path)}, skipping.")

# Process each .fasta.gz file if a corresponding .tsv file exists
for fasta_file_path in glob.glob(os.path.join(fasta_dir_path, '*.fasta.gz')):
    base_name = os.path.basename(fasta_file_path).replace('.fasta.gz', '')
    gene_info_path = os.path.join(tsv_dir_path, f"{base_name}_{gene_name}_gene.tsv")
    output_fasta_path = os.path.join(modified_dir_path, f"{base_name}_modified.fasta.gz")

    if os.path.exists(gene_info_path):
        gene_info = load_gene_info(gene_info_path)
        modify_sequences(fasta_file_path, gene_info, output_fasta_path)
    else:
        print(f"No gene info file found for {fasta_file_path}, skipping modification.")