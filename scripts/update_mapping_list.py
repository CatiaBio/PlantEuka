import os
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Update mapping list.')
parser.add_argument('--mapping_file', required=True, help='Path to mapping file')
parser.add_argument('--input_folder', required=True, help='Path to input folder')
parser.add_argument('--output_file', required=True, help='Path to the output directory')
args = parser.parse_args()

# Use argparse to parse command line arguments
mapping_file_path = args.mapping_file
input_folder_path = args.input_folder
output_dir_path = args.output_file

# List all files in the specified folder
files_in_folder = os.listdir(input_folder_path)

# Prepare a list to hold lines for IDs not found in the folder
lines_not_in_folder = []

with open(mapping_file_path, 'r') as file:
    for line in file:
        ncbi_id = line.strip().split('\t')[0]
        # Check for the presence of the NCBI ID in the folder, considering the file extension .fasta.gz
        file_name = f"{ncbi_id}.fasta.gz"  # Adjusted to match the file naming convention
        if file_name not in files_in_folder:
            lines_not_in_folder.append(line.strip())  # Keep the whole line including the taxid
        # else:
        #     print(f"Found {ncbi_id} in the folder. Removing from list.")


# lines_not_in_folder now contains lines for IDs not found in the folder, including their taxids
#print("Lines for NCBI IDs not found in the folder (including taxids):")
# for line in lines_not_in_folder:
#     print(line)

with open(output_dir_path, 'w') as file:
    last_index = len(lines_not_in_folder) - 1  # Get the index of the last item
    for index, line in enumerate(lines_not_in_folder):
        if index != last_index:
            file.write(f"{line}\n")  # Write with a newline for all but the last line
        else:
            file.write(line)  # Write without adding a newline for the last line