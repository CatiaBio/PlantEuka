#!/bin/bash

# Define the base directory containing the taxonomy folders
base_dir="genomes/chloroplast/sorted"

# Define the directory where pair lists will be saved
output_dir="other/lists"

# Ensure the output directory exists
mkdir -p "$output_dir"
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
        local output_file="${output_dir}/${taxonomy_level}_${folder_name}_pairs.txt"
        local list_file_path="${folder}/list_of_fasta_files.txt"

        # Generate a list of all fasta.gz files
        find "$folder" -name '*.fasta.gz' > "$list_file_path"
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

        > "$output_file"  # Clear or create the output file
        for (( i=0; i<${#files[@]}; i++ )); do
            for (( j=i+1; j<${#files[@]}; j++ )); do
                echo "${files[i]} ${files[j]}" >> "$output_file"
            done
        done

        # Clean up
        rm "$list_file_path"  # Optionally remove the temporary list file
        echo "Generated pairs list saved to $output_file"
    done
}

# Process each taxonomy level
for level in genus family order; do
    process_folders "$level"
done