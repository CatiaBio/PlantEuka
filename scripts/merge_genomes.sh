#!/bin/bash

# Check if the correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <base_directory> <merged_directory>"
    exit 1
fi

# Assign command-line arguments to variables
base_dir="$1"
merged_dir="$2"

# Ensure the merged directory and its subdirectory 'acc_numb' exist
mkdir -p "$merged_dir"
mkdir -p "${merged_dir}/acc_numb"  # Create 'acc_numb' subdirectory

# Function to process folders and merge FASTA files
process_folder() {
    local folder="$1"
    local category=$(basename "$(dirname "$folder")")
    local sub_dir_name=$(basename "$folder")

    # Skip processing for 'original_before_clean' directory
    if [[ "$sub_dir_name" == "original_before_clean" ]]; then
        return
    fi

    # Output file names
    local output_file="${merged_dir}/${category}_${sub_dir_name}.fasta.gz"
    # Modify tsv_file path to save inside the 'acc_numb' folder
    local tsv_file="${merged_dir}/acc_numb/${category}_${sub_dir_name}.tsv"  # Adjusted path

    # Find and merge FASTA files
    find "$folder" -type f -name "*.fasta.gz" > /dev/null
    if [ $? -eq 0 ]; then
        find "$folder" -type f -name "*.fasta.gz" -print0 | xargs -0 zcat | gzip > "$output_file"
        # Extract acc_numb after confirming file creation
        zcat "$output_file" | grep "^>" | sed 's/>//g' | cut -d ' ' -f 1 > "$tsv_file"
    fi
}

# Categories to process
categories=("genus")
#categories=("genus" "family" "order")

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
