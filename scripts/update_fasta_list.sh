#!/bin/bash

# Define the base directory where the folders are located
BASE_DIR="data/mitochondrion/organized"

# Define the file containing your list
LIST_FILE="data/mapping_id_species_mt.txt"

# Temporary file to store updated list
UPDATED_LIST="data/updated_list_mt.txt"

# Clear the updated list file to start with an empty file
> "$UPDATED_LIST"

# Read each line from the input list
while IFS=$'\t' read -r fasta_id species; do
    # Initialize a flag to indicate if the file was found
    file_found=false

    # Construct the file path pattern for recursive search within genus, family, and order directories
    FILE_PATTERN="${BASE_DIR}/*/*/${fasta_id}.fasta.gz"
    
    # Check if the file exists in any of the specified paths
    if ls $FILE_PATTERN 1> /dev/null 2>&1; then
        file_found=true
        # If the file is found, append the line to the updated list
        echo -e "${fasta_id}\t${species}" >> "$UPDATED_LIST"
    else
        echo "File not found for ${fasta_id}, ${species}"
    fi
done < "$LIST_FILE"

# Optionally, replace the original list with the updated list
# mv "$UPDATED_LIST" "$LIST_FILE"

echo "Updated list saved to $UPDATED_LIST"
