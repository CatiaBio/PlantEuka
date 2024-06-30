import csv
import glob
import os
import sys

# Check for command-line arguments for gene name
if len(sys.argv) < 2:
    print("Usage: python script_name.py <gene_name>")
    sys.exit(1)  # Exit the script if gene name is not provided

gene_name = sys.argv[1]  # Get the gene name from the command-line argument

# Specify the file path for the gene data dynamically based on the gene name
test_gene_file_path = f'genes/chloroplast/{gene_name}/{gene_name}_gene_info.tsv'

# Directory containing the .tsv files to process
input_dir_path = 'genomes/chloroplast/merged/acc_numb'

# Define the base directory for results dynamically based on the gene name
results_base_dir = f'genomes/chloroplast/merged/{gene_name}'

# Define the directory for missing files dynamically based on the gene name
missing_files_dir = f'genomes/chloroplast/merged/{gene_name}/missing'

# Ensure the directories exist
os.makedirs(results_base_dir, exist_ok=True)
os.makedirs(missing_files_dir, exist_ok=True)

# List all .tsv files in the directory
tsv_files = glob.glob(os.path.join(input_dir_path, '*.tsv'))

# Function to process each .tsv file
def process_file(accession_file_path, test_clpP_file_path):
    base_name = os.path.basename(accession_file_path).replace('.tsv', '')
    output_file_path = os.path.join(results_base_dir, f"{base_name}_{gene_name}_gene.tsv")
    missing_file_path = os.path.join(missing_files_dir, f"{base_name}_missing.tsv")
    
    # Read accession numbers from the tsv file
    accession_numbers = set()
    with open(accession_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            accession_numbers.add(row[0].split(' ')[0])

    # Prepare to collect matching entries
    matching_entries = []
    found_accessions = set()

    # Read and process the test clpP file for matches
    with open(test_clpP_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            acc_version, location_info = row[0].split(' ', 1)
            strand = location_info.split(']')[-1]  # Extract the strand information
            if acc_version in accession_numbers:
                start_end = location_info.split(']')[0].strip('[')
                start, end = start_end.split(':')
                matching_entries.append({'acc_version': acc_version, 'start_pos': int(start) + 1, 'end_pos': int(end), 'orientation': strand})
                found_accessions.add(acc_version)

    # Identify missing accession numbers
    missing_accessions = accession_numbers - found_accessions

    # Write matched entries to the output file, if all accession numbers are found
    if not missing_accessions:
        with open(output_file_path, 'w', newline='') as output_file:
            fieldnames = ['acc_version', 'start_pos', 'end_pos', 'orientation']
            writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for entry in matching_entries:
                writer.writerow(entry)
    else:
        # Write missing accession numbers to a file
        with open(missing_file_path, 'w', newline='') as missing_file:
            writer = csv.writer(missing_file, delimiter='\t')
            for acc in missing_accessions:
                writer.writerow([acc])
        print(f"Created missing accession numbers file: {missing_file_path}")

# Process each .tsv file in the directory
for tsv_file in tsv_files:
    if f'_{gene_name}_gene.tsv' not in tsv_file:  # Dynamically skip files based on gene name
        process_file(tsv_file, test_gene_file_path)
