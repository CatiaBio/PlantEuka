#!/bin/bash

# Check if the correct number of command-line arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <organelle>"
    exit 1
fi

# Extract command-line arguments
organelle="$1"

# Define the base directory containing the taxonomy folders
base_dir="$(pwd)/${organelle}/genomes/sorted"

# Define the directories for the outputs and logs
output_dir="$(pwd)/${organelle}/other/pairwise"
results_dir="$(pwd)/${organelle}/results/pairwise"

# Ensure the output and results directories exist
mkdir -p "$output_dir" "$results_dir"
if [ ! -d "$output_dir" ]; then
    echo "Failed to create or access the output directory: $output_dir"
    exit 1
fi

# Function to process folders of a given taxonomic level
process_folders() {
    local taxonomy_level="$1"  # e.g., genus, family, or order
    local taxonomy_dir="${base_dir}/${taxonomy_level}"

    echo "Processing $taxonomy_level directories in $taxonomy_dir"

    # Iterate over the folders at the given taxonomy level
    for folder in "$taxonomy_dir"/*; do
        local folder_name=$(basename "$folder")

        # Skip 'original_before_clean' directory and verify directory existence
        if [[ "$folder_name" == "original_before_clean" ]] || [ ! -d "$folder" ]; then
            echo "Skipping $folder"
            continue
        fi

        echo "Processing folder: $folder"
        local list_file_path="${folder}/list_of_fasta_files.txt"

        # Generate a list of all fasta.gz files with full paths
        find "$folder" -type f -name '*.fasta.gz' > "$list_file_path"
        if [ ! -s "$list_file_path" ]; then
            echo "No FASTA files found in $folder"
            continue
        fi

        # Read the list file and generate pairs
        exec 3< "$list_file_path"
        files=()
        while IFS= read -r -u 3 file; do
            files+=("$file")
        done

        for (( i=0; i<${#files[@]}; i++ )); do
            for (( j=i+1; j<${#files[@]}; j++ )); do
                pair1="${files[i]}"
                pair2="${files[j]}"
                path_entry=(${pair1//\// })
                taxonomy_genus=$(echo "${path_entry[6]}")
                taxonomy_species=$(echo "${path_entry[7]}")
                pair1_name=$(basename "${pair1}" .fasta.gz)
                pair2_name=$(basename "${pair2}" .fasta.gz)
                output_name="${results_dir}/${taxonomy_level}_${folder_name}_${pair1_name}_${pair2_name}.stretcher"
                echo "$pair1 $pair2 $output_name" >> "${output_dir}/${taxonomy_level}_${folder_name}_pairs.txt"
            done
        done

        # Optionally remove the temporary list file
        rm "$list_file_path"
        echo "Generated pairs list saved to ${output_dir}/${taxonomy_level}_${folder_name}_pairs.txt"
    done
}

# Process each taxonomy level
for level in genus family order; do
    process_folders "$level"
done

# Concatenate all pair list files into one file
cat ${output_dir}/*_pairs.txt > ${output_dir}/all_pairs_parallel.txt
echo "All pair lists concatenated into ${output_dir}/all_pairs_parallel.txt"
