import os

def list_files_with_results(base_folder, input_subfolder, results_subfolder, output_file_subfolder, output_filename):
    # Create the full paths based on the base folder
    input_folder = os.path.join(base_folder, input_subfolder)
    results_folder = os.path.join(base_folder, results_subfolder)
    output_file = os.path.join(base_folder, output_file_subfolder, output_filename)
    
    # Create a list to hold the final output
    file_list = []

    # Iterate over all files in the input folder
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            # Construct the full file path
            file_path = os.path.join(root, file)
            
            # Determine the result file path with .aln extension
            file_name_without_fasta_gz = file.replace('.fasta.gz', '')  # Remove .fasta.gz extension
            result_file_name = f"{file_name_without_fasta_gz}.aln"  # Add .aln extension
            result_file_path = os.path.join(results_folder, result_file_name)
            
            # Add the file paths to the list
            file_list.append(f"{file_path} {result_file_path}")

    # Save the list to a file
    with open(output_file, 'w') as f:
        for file_pair in file_list:
            f.write(file_pair + "\n")

# Determine the base folder dynamically (assumes the PlantEuka folder is in the user's home directory)
home_dir = os.path.expanduser("~")
base_folder = os.path.join(home_dir, 'Projects', 'PlantEuka')

# Define the subfolders and output file name
input_subfolder = 'chloroplast/genomes/merged'
results_subfolder = 'chloroplast/results/msa'
output_file_subfolder = 'chloroplast/other'
output_filename = 'msa_paths.txt'

# Generate the list of files with result paths and save to the output file
list_files_with_results(base_folder, input_subfolder, results_subfolder, output_file_subfolder, output_filename)
