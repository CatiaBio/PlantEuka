#!/bin/bash

# Define the size of the dummy file (1 KB)
dummy_size=1024  # Size in bytes

# Find all files ending with .stretcher.gz and create dummy files with .stretcher extension
find . -type f -name "*.stretcher.gz" | while read -r file; do
    # Get the base name without the .gz extension
    base_name="${file%.gz}"
    # Create the dummy file with .stretcher extension
    dummy_file="${base_name%.gz}"
    # Create a dummy file with the specified size
    head -c $dummy_size /dev/urandom > "$dummy_file"
    echo "Created dummy file: $dummy_file"
done

echo "All dummy files created."