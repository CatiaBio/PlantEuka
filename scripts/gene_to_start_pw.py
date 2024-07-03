#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
import csv
import gzip
import glob
import os
import sys

fasta_dir_path = 'chloroplast/genomes/original'
tsv_file_path = 'chloroplast/other/accession_with_gene_info.txt'
merged_dir_path = 'chloroplast/genomes/merged'
modified_dir_path = 'genomes/chloroplast/merged/modified'

os.makedirs(merged_dir_path, exist_ok=True)  # Ensure the merged directory exists
os.makedirs(modified_dir_path, exist_ok=True)  # Ensure the modified directory exists

def load_gene_info(tsv_file_path):
    """Load start and end positions and orientation for genes from a TSV file."""
    gene_info = {}
    current_key = None
    with open(tsv_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith("order_") or row[0].startswith("family_"):
                current_key = row[0]
                gene_info[current_key] = []
            else:
                accession, _, start_pos, end_pos, orientation = row
                gene_info[current_key].append((accession, int(start_pos), int(end_pos), orientation))
    return gene_info

def modify_sequences(fasta_file_path, gene_info, output_fasta_path):
    """Modify sequences by applying reverse complement if necessary, then rotate to start at gene's position."""
    with gzip.open(fasta_file_path, 'rt') as file, gzip.open(output_fasta_path, 'wt') as output_file:
        for record in SeqIO.parse(file, 'fasta'):
            acc_version = record.id.split('.')[0]  # Use the accession number without version
            for info in gene_info:
                if acc_version == info[0]:
                    start_pos, end_pos, orientation = info[1], info[2], info[3]

                    # Apply reverse complement first if orientation is negative
                    if orientation == 'minus':
                        record.seq = record.seq.reverse_complement()
                        # Adjust start and end positions for the new sequence orientation
                        start_pos, end_pos = len(record.seq) - end_pos + 1, len(record.seq) - start_pos + 1

                    # Now, rotate the sequence to start at the gene's start position
                    rotated_seq = record.seq[start_pos-1:] + record.seq[:start_pos-1]
                    record.seq = rotated_seq

                    SeqIO.write([record], output_file, 'fasta')
                    break
            else:
                print(f"Gene info not found for {acc_version} in file {os.path.basename(fasta_file_path)}, skipping.")

def merge_and_modify_sequences(order, gene_info_list):
    """Merge and modify sequences for a given order."""
    merged_fasta_path = os.path.join(merged_dir_path, f"{order}.fasta.gz")
    modified_fasta_path = os.path.join(modified_dir_path, f"{order}_modified.fasta.gz")

    with gzip.open(merged_fasta_path, 'wt') as merged_file:
        for gene_info in gene_info_list:
            accession = gene_info[0]
            fasta_file_path = os.path.join(fasta_dir_path, f"{accession}.fasta.gz")
            
            if os.path.exists(fasta_file_path):
                with gzip.open(fasta_file_path, 'rt') as fasta_file:
                    for record in SeqIO.parse(fasta_file, 'fasta'):
                        SeqIO.write([record], merged_file, 'fasta')
            else:
                print(f"No fasta file found for {accession}, skipping merge.")

    # Now modify the merged file
    if os.path.exists(merged_fasta_path):
        modify_sequences(merged_fasta_path, gene_info_list, modified_fasta_path)
    else:
        print(f"No merged fasta file found for {order}, skipping modification.")

def main(tsv_file_path):
    """Main function to handle the merging and modification of sequences."""
    gene_info = load_gene_info(tsv_file_path)
    
    for order, gene_info_list in gene_info.items():
        merge_and_modify_sequences(order, gene_info_list)

if __name__ == "__main__":
    main(tsv_file_path)
