#!/bin/bash

# Description:
# This script automates the process of "cleaning" FASTA files located in specified directories.
# It searches for FASTA files (*.fasta.gz), replaces specified nucleotide ambiguity codes with 'N',
# and saves the cleaned sequences in new files. The script logs its actions to a user-specified log file.
# Usage:
# ./script_name.sh <base_directory> <log_directory> <log_file_name>
# Example:
# ./scripts/clean_genomes.sh genomes/mitochondrion/sorted logs/genome_cleanup_mt.log

# Check if the correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <base_directory> <log_file_name>"
    exit 1
fi

# Assign command-line arguments to variables
base_dir="$1"
log_file_name="$2"

# Start new log file or overwrite if already exists
echo "Starting cleanup process at $(date)" > "$log_file_name"

# Categories to process
categories=("genus" "family" "order")

# Loop through each category
for category in "${categories[@]}"; do
    # Construct the path to the current category
    category_path="${base_dir}/${category}"

    if [ -d "$category_path" ]; then
        echo "Processing category: $category" >> "$log_file_name"
        # Find all FASTA.gz files in the current category
        find "$category_path" -type f -name "*.fasta.gz" | while read fasta_file; do
            # Extract the base name without the .fa.gz extension and directory
            base_name=$(basename "$fasta_file" .fa.gz)
            
            # Construct the path for the cleaned FASTA file in the same directory
            cleaned_fasta="${fasta_file%.*}_cleaned.fa.gz"
            
            # Read and clean the FASTA file's content, replacing specific codes with 'N'
            original_content=$(zcat "$fasta_file")
            cleaned_content=$(echo "$original_content" | sed '/^>/!s/[RYSWKSMDVHBX]/N/g')

            # Compare the original and cleaned contents
            if [ "$original_content" != "$cleaned_content" ]; then
                # If different, save the cleaned content and log the action
                echo "$cleaned_content" | gzip > "$cleaned_fasta"
                echo "Cleaned FASTA file created: $cleaned_fasta" >> "$log_file_name"
            else
                # If no changes were made, log that cleaning was not necessary
                echo "No cleaning needed for $fasta_file" >> "$log_file_name"
            fi
        done
    else
        # Log if the directory for a category does not exist
        echo "Directory does not exist: $category_path" >> "$log_file_name"
    fi
done

# Log completion of the process
echo "All analyses completed at $(date)" >> "$log_file_name"
