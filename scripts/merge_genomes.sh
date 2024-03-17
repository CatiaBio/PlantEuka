#!/bin/bash

# Check if the correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <base_directory> <merged_directory>"
    exit 1
fi

# Assign command-line arguments to variables
base_dir="$1"
merged_dir="$2"

# Ensure the merged directory exists
mkdir -p "$merged_dir"

# Categories to process
categories=("genus" "family" "order")

# Loop through each category
for category in "${categories[@]}"; do
    # Path to the current category
    category_path="${base_dir}/${category}"

    # Check if the category directory exists
    if [ -d "$category_path" ]; then
        # Find all subdirectories in the current category
        for sub_dir in "$category_path"/*; do
            if [ -d "$sub_dir" ]; then
                # Extract the name of the subdirectory
                sub_dir_name=$(basename "$sub_dir")

                # Define the output file names
                output_file="${merged_dir}/${category}_${sub_dir_name}.fa.gz"
                tsv_file="${merged_dir}/${category}_${sub_dir_name}.tsv"

                # Initialize the TSV file
                > "$tsv_file"

                # Prepare to merge, prioritizing '_cleaned' files
                for fasta_file in "$sub_dir"/*_cleaned.fa.gz; do
                    # Extract ID (file name without path and extension)
                    id=$(basename "$fasta_file" "_cleaned.fa.gz")
                    echo -e "$id" >> "$tsv_file"
                done

                # If no '_cleaned' files found, fallback to original '.fa.gz' files
                if [ ! -s "$tsv_file" ]; then
                    for fasta_file in "$sub_dir"/*.fa.gz; do
                        # Extract ID (file name without path and extension)
                        id=$(basename "$fasta_file" .fa.gz)
                        echo -e "$id" >> "$tsv_file"
                    done
                fi

                # Merge and compress selected FASTA files into one
                zcat "$sub_dir"/*_cleaned.fa.gz "$sub_dir"/*.fa.gz 2>/dev/null | gzip > "$output_file"

                echo "Merged file created: $output_file"
                echo "TSV file created: $tsv_file"
            fi
        done
    else
        echo "Directory does not exist: $category_path"
    fi
done
