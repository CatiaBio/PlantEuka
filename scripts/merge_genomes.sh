#!/bin/bash

# Description: 
# This script merges FASTA files within subdirectories of a base directory.

# Check if the correct number of arguments are passed
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <organelle>"
    exit 1
fi

# Assign command-line argument to a variable
organelle="$1"

# Define base and merged directories based on the organelle
base_dir="${organelle}/genomes/sorted"
merged_dir="${organelle}/genomes/merged"

# Ensure the merged directory exists
mkdir -p "$merged_dir"

# Function to process folders and merge FASTA files
process_folder() {
    local folder="$1"
    local category=$(basename "$(dirname "$folder")")
    local sub_dir_name=$(basename "$folder")

    # Skip processing for 'original_before_clean' directory
    if [[ "$sub_dir_name" == "original_before_clean" ]]; then
        return
    fi

    # Output file name
    local output_file="${merged_dir}/${category}_${sub_dir_name}.fasta.gz"

    # Find and merge FASTA files
    find "$folder" -type f -name "*.fasta.gz" > /dev/null
    if [ $? -eq 0 ]; then
        find "$folder" -type f -name "*.fasta.gz" -print0 | xargs -0 zcat | gzip > "$output_file"
    fi
}

# Categories to process
categories=("genus" "family" "order")

# Loop through each category
for category in "${categories[@]}"; do
    # Path to the current category
    category_path="${base_dir}/${category}"

    # Check if the category directory exists
    if [ -d "$category_path" ]; then
        echo "Merging $category"
        # Find all subdirectories within the category directory
        for dir in "$category_path"/*/; do
            if [ -d "$dir" ] && [[ "$(basename "$dir")" != "original_before_clean" ]]; then
                process_folder "$dir"
            fi
        done
        echo "Done"
    fi
done