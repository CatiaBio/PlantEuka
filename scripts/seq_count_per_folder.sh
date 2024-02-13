#!/bin/bash

# Specify output files
output_file1="../data/mitochondrion/folder_file_counts_mt.txt"
output_file2="../data/chloroplast/folder_file_counts_cp.txt"

# Specify the parent folders
parent_folder1="../data/mitochondrion/"
parent_folder2="../data/chloroplast/"

# Function to process folders
process_folder() {
    local parent_folder=$1
    local output_file=$2

    # Iterate through each subfolder and count the files
    for folder in "$parent_folder"/*/; do
        folder_name=$(basename "$folder")

        # Remove "_mt" or "_cp" suffix from the folder name
        folder_name_without_suffix=${folder_name%_mt}
        folder_name_without_suffix=${folder_name_without_suffix%_cp}

        file_count=$(find "$folder" -type f | wc -l)
        echo "$folder_name_without_suffix: $file_count" >> "$output_file"
    done
}

# Process mitochondrion folder
process_folder "$parent_folder1" "$output_file1"

# Process chloroplast folder
process_folder "$parent_folder2" "$output_file2"
