from Bio.Seq import Seq
from Bio import SeqIO
import csv
import gzip
import glob
import os
import sys

# Check for command-line arguments for gene name
if len(sys.argv) < 2:
    print("Usage: python script_name.py <gene_name>")
    sys.exit(1)

gene_name = sys.argv[1]  # Gene name from command-line argument

fasta_dir_path = 'genomes/chloroplast/merged'
tsv_dir_path = f'genomes/chloroplast/merged/{gene_name}'
modified_dir_path = f'genomes/chloroplast/merged/modified_{gene_name}'

os.makedirs(modified_dir_path, exist_ok=True)  # Ensure the modified directory exists

def load_gene_info(gene_info_path):
    """Load start and end positions and orientation for genes from a TSV file."""
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
    """Modify sequences by applying reverse complement if necessary, then rotate to start at gene's position."""
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
                # For a reverse complemented sequence, the start and end positions have been adjusted above
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
