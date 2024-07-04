#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
import csv
import gzip
import os
import re

fasta_dir_path = 'chloroplast/genomes/original'
tsv_file_path = 'chloroplast/other/accession_with_gene_info_original.txt'
modified_dir_path = 'chloroplast/genomes/modified_with_gene'
merged_dir_path = 'chloroplast/genomes/merged'

os.makedirs(merged_dir_path, exist_ok=True)  # Ensure the merged directory exists
os.makedirs(modified_dir_path, exist_ok=True)  # Ensure the modified directory exists

def load_gene_info(tsv_file_path):
    """Load start and end positions and orientation for genes from a TSV file."""
    gene_info = []
    current_group = None
    accession_set = set()
    with open(tsv_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith(('family', 'order', 'genus')):
                current_group = row[0]
                print(f"Starting with {current_group}")
            else:
                if current_group:
                    accession = row[0]
                    gene = row[1]
                    start_pos = int(row[2])
                    end_pos = int(row[3])
                    orientation = row[4]
                    gene_info.append((current_group, accession, gene, start_pos, end_pos, orientation))
                    accession_set.add(accession)
                    print(f"Processing accession {accession}")
    return gene_info, accession_set

def modify_sequences(fasta_file_path, gene_info, output_fasta_path):
    """Modify sequences by applying reverse complement if necessary, then rotate to start at gene's position."""
    accession_pattern = re.compile(r'^[A-Za-z0-9_.]+')
    with gzip.open(fasta_file_path, 'rt') as file, gzip.open(output_fasta_path, 'wt') as output_file:
        for record in SeqIO.parse(file, 'fasta'):
            acc_version_match = accession_pattern.match(record.id)
            if acc_version_match:
                acc_version = acc_version_match.group(0)
                for info in gene_info:
                    if acc_version == info[1]:
                        start_pos, end_pos, orientation = info[3], info[4], info[5]

                        # Apply reverse complement first if orientation is minus
                        if orientation == 'minus':
                            record.seq = record.seq.reverse_complement()
                            # Adjust start and end positions for the new sequence orientation
                            start_pos, end_pos = len(record.seq) - end_pos, len(record.seq) - start_pos

                        # Now, rotate the sequence to start at the gene's start position
                        rotated_seq = record.seq[start_pos:end_pos] + record.seq[:start_pos] + record.seq[end_pos:]
                        record.seq = rotated_seq

                        SeqIO.write([record], output_file, 'fasta')
                        break

def modify_all_sequences(fasta_dir_path, gene_info, modified_dir_path):
    """Modify all sequences in the fasta directory according to the gene info."""
    for gene_info_item in gene_info:
        accession = gene_info_item[1]
        fasta_file_path = os.path.join(fasta_dir_path, f"{accession}.fasta.gz")
        output_fasta_path = os.path.join(modified_dir_path, f"{accession}_modified.fasta.gz")
        
        if os.path.exists(fasta_file_path):
            print(f"Modifying {accession}")
            modify_sequences(fasta_file_path, gene_info, output_fasta_path)

def merge_sequences(group_name, gene_info_list, modified_dir_path, merged_dir_path):
    """Merge modified sequences for a given family, order, or genus."""
    merged_fasta_path = os.path.join(merged_dir_path, f"{group_name.replace(' ', '_')}.fasta.gz")
    
    with gzip.open(merged_fasta_path, 'wt') as merged_file:
        for gene_info in gene_info_list:
            accession = gene_info[1]
            modified_fasta_path = os.path.join(modified_dir_path, f"{accession}_modified.fasta.gz")
            
            if os.path.exists(modified_fasta_path):
                with gzip.open(modified_fasta_path, 'rt') as fasta_file:
                    for record in SeqIO.parse(fasta_file, 'fasta'):
                        SeqIO.write([record], merged_file, 'fasta')
    print(f"Done with {group_name}")

def main(tsv_file_path, fasta_dir_path, modified_dir_path, merged_dir_path):
    """Main function to handle the modification and merging of sequences."""
    gene_info, accession_set = load_gene_info(tsv_file_path)
    
    # Step 1: Modify all sequences and save them in the modified directory
    #modify_all_sequences(fasta_dir_path, gene_info, modified_dir_path)
    
    # Step 2: Group gene info by family, order, or genus and merge the modified sequences
    grouped_info = {}
    for info in gene_info:
        group_name = info[0]
        if group_name not in grouped_info:
            grouped_info[group_name] = []
        grouped_info[group_name].append(info)
    
    for group_name, gene_info_list in grouped_info.items():
        merge_sequences(group_name, gene_info_list, modified_dir_path, merged_dir_path)

if __name__ == "__main__":
    main(tsv_file_path, fasta_dir_path, modified_dir_path, merged_dir_path)