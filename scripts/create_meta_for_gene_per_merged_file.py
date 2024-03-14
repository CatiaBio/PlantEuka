import os
import csv

# Path to the folder containing your .tsv files
folder_path = '../data/chloroplast/merged'  # Replace with your folder path

# Path to your metadata file
metadata_file_path = '../other_info/rubisco_gene_meta.txt'  # Replace with the path to your metadata file

# Function to perform matching and create new file
def match_and_create_file(tsv_file, metadata_file_path):
    # Read the accessions into a set for fast lookup
    accessions = set()
    with open(tsv_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            accessions.add(row[0])  # Assuming the accession is the first element

    # New filename for matched records
    output_file_path = tsv_file.replace('.tsv', '_gene_rbcL.tsv')

    # Counter for missing matches
    missing_matches = 0

    # Read metadata file and write out the matching records to the new output file
    with open(metadata_file_path, 'r') as metadata_file, open(output_file_path, 'w', newline='') as output_file:
        metadata_reader = csv.reader(metadata_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')

        for row in metadata_reader:
            # Check if the 'genomic_nucleotide_accession.version' column matches any accession
            if row[11] in accessions:  # Adjust index 8 if necessary
                writer.writerow(row)
            else:
                missing_matches += 1

# Print message about matched records
    print(f"Matching records from {tsv_file} have been saved to: {output_file_path}")
    if missing_matches > 0:
        print(f"Warning: There are {missing_matches} accessions from {tsv_file} that were not found in the metadata.")

# Iterate over .tsv files in the specified folder and perform the match operation
for filename in os.listdir(folder_path):
    if filename.endswith('.tsv'):
        file_path = os.path.join(folder_path, filename)
        match_and_create_file(file_path, metadata_file_path)
