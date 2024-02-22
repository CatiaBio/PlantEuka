import os
import gzip

# Define the base directory where the organized FASTA files are located
base_dir = 'data/mitochondrion/organized'  # Update this path to your specific directory

# Initialize a list to hold all the paths of FASTA files
fasta_files = []

# Walk through the directory structure starting from base_dir
for root, dirs, files in os.walk(base_dir):
    for file in files:
        # Check if the file is a .fasta.gz file
        if file.endswith('.fasta.gz'):
            # Add the full path of the file to the list
            fasta_files.append(os.path.join(root, file))

# Path to the output file that will contain all FASTA sequences combined
output_file_path = 'data/combined_sequences_mt.fa.gz'  # Update this path as needed

# Combine all FASTA files into one
with gzip.open(output_file_path, 'wb') as outfile:
    for fasta_file in fasta_files:
        with gzip.open(fasta_file, 'rb') as infile:
            outfile.write(infile.read())

print(f"Combined FASTA file created at: {output_file_path}")
