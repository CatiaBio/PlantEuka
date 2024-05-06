#!/bin/bash

# Define the base directory containing the taxonomy folders
base_dir="genomes/chloroplast/sorted"  # Change this to your base directory

# Function to find and unzip all .gz files in a given directory, excluding 'original_before_clean'
unzip_files() {
    local directory=$1
    echo "Searching in $directory for .gz files to unzip, excluding any 'original_before_clean' directories..."

    # Find all .gz files excluding those in 'original_before_clean' directories
    find "$directory" -type f -name '*.gz' -not -path '*/original_before_clean/*' | while read gz_file; do
        local output_file="${gz_file%.gz}"  # Remove the .gz extension for the output file name
        
        # Check if the output file already exists
        if [ ! -f "$output_file" ]; then
            # Unzip the file while keeping the original
            gzip -k -d "$gz_file"
            echo "Unzipped $gz_file"
        else
            echo "Output file $output_file already exists, skipping."
        fi
    done

    echo "Unzipping complete in $directory."
}

# Unzip files in genus, family, and order directories
for taxonomy in genus family order; do
    taxonomy_dir="${base_dir}/${taxonomy}"
    if [ -d "$taxonomy_dir" ]; then
        unzip_files "$taxonomy_dir"
    else
        echo "Directory $taxonomy_dir does not exist."
    fi
done