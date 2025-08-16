import os

# For Snakemake compatibility
input_folder = snakemake.input.merged_dir
output_file = snakemake.output.msa_paths
organelle = snakemake.params.organelle

# Create results folder path
results_folder = f"{organelle}/results/msa"

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)
os.makedirs(results_folder, exist_ok=True)

# Create a list to hold the final output
file_list = []

# Iterate over all files in the input folder
for root, dirs, files in os.walk(input_folder):
    for file in files:
        if file.endswith('.fasta.gz'):
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

print(f"Generated MSA paths for {len(file_list)} files")
