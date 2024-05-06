from Bio.Seq import Seq
from Bio import SeqIO
import csv
import gzip
import glob
import os
import sys
from collections import Counter
from datetime import datetime

def load_gene_info(gene_info_path):
    """Load gene information from TSV file."""
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

def clean_sequence(sequence):
    """Clean sequence by replacing non-nucleotide characters with 'N'."""
    replaced_details = []  
    cleaned_sequence = []
    for index, nucleotide in enumerate(sequence.upper()):
        if nucleotide not in 'ACGTN':
            cleaned_sequence.append('N')
            replaced_details.append(nucleotide)  
        else:
            cleaned_sequence.append(nucleotide)
    return ''.join(cleaned_sequence), Counter(replaced_details)

def modify_and_clean_sequences(fasta_file_path, gene_info, output_fasta_path):
    """Modify and clean sequences, return log entries for cleaning and reverse complementation."""
    clean_log_entries = []
    reverse_complement_log_entries = []
    modified_sequences = []
    with gzip.open(fasta_file_path, 'rt') as file:
        for record in SeqIO.parse(file, 'fasta'):
            acc_version = record.id.split(' ')[0]
            if acc_version in gene_info:
                start_pos, end_pos, orientation = gene_info[acc_version]
                if orientation == '(-)':
                    record.seq = record.seq.reverse_complement()
                    reverse_complement_log_entries.append(f"{record.id}\n")
                
                rotated_seq = record.seq[start_pos-1:] + record.seq[:start_pos-1]
                cleaned_seq, replacements = clean_sequence(str(rotated_seq))
                record.seq = Seq(cleaned_seq)
                
                if replacements:
                    replaced_info = ', '.join(f"{char}({count})" for char, count in replacements.items())
                    clean_log_entries.append(f"{record.id}\t{sum(replacements.values())}\t{replaced_info}\n")
                    
                modified_sequences.append(record)

    with gzip.open(output_fasta_path, 'wt') as output_file:
        SeqIO.write(modified_sequences, output_file, 'fasta')
    
    return clean_log_entries, reverse_complement_log_entries

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <gene_name>")
        sys.exit(1)

    gene_name = sys.argv[1]
    fasta_dir_path = 'genomes/chloroplast/merged'
    tsv_dir_path = f'genomes/chloroplast/merged/{gene_name}'
    modified_dir_path = f'genomes/chloroplast/merged/modified_{gene_name}'
    os.makedirs(modified_dir_path, exist_ok=True)

    clean_log_path = os.path.join(modified_dir_path, "clean.log")
    reverse_complement_log_path = os.path.join(modified_dir_path, "reverse_complement.log")

    all_clean_log_entries = []
    all_reverse_complement_log_entries = []

    for fasta_file_path in glob.glob(os.path.join(fasta_dir_path, '*.fasta.gz')):
        base_name = os.path.basename(fasta_file_path).replace('.fasta.gz', '')
        gene_info_path = os.path.join(tsv_dir_path, f"{base_name}_{gene_name}_gene.tsv")

        if os.path.exists(gene_info_path):
            gene_info = load_gene_info(gene_info_path)
            output_fasta_path = os.path.join(modified_dir_path, f"{base_name}_modified.fasta.gz")
            clean_log_entries, reverse_complement_log_entries = modify_and_clean_sequences(fasta_file_path, gene_info, output_fasta_path)
            all_clean_log_entries.extend(clean_log_entries)
            all_reverse_complement_log_entries.extend(reverse_complement_log_entries)
        else:
            print(f"Gene info file not found for {base_name}, skipping.")

    # Writing logs at the end
    with open(clean_log_path, 'w') as clean_log:
        clean_log.write("ID\tNumber_Replacements\tReplaced_Characters\n")
        for entry in all_clean_log_entries:
            clean_log.write(entry)

    with open(reverse_complement_log_path, 'w') as reverse_complement_log:
        reverse_complement_log.write("ID\n")
        for entry in all_reverse_complement_log_entries:
            reverse_complement_log.write(entry)