#!/bin/bash

# Path to the list of genome IDs
LIST_PATH="/home/projects/MAAG/PlantEuka/other/afterpw_genomes_to_remove_cp.tsv"

# Base directory containing the subfolders order, genus, and family
BASE_DIR="/home/projects/MAAG/PlantEuka/genomes/chloroplast/sorted_after_pw"

# Read each line from the list
while read -r genome_id; do
  # Find and remove the corresponding .fasta.gz files in the directory tree
  find "$BASE_DIR" -type f -name "${genome_id}.fasta.gz" -exec rm {} +
done < "$LIST_PATH"

echo "Deletion complete."
