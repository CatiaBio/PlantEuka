#!/bin/bash

# Function to extract taxid using zgrep and awk
extract_taxid() {
  accession_version=$1
  temp_db_file=$2
  output_file=$3
  result=$(grep -P "^${accession_version}\t" "$temp_db_file" | awk -v acc="$accession_version" '$1 == acc {print $1 "\t" $2}')
  if [ -n "$result" ]; then
    echo -e "$result" >> "$output_file"
  else
    echo -e "${accession_version}\tNot_Found" >> "$output_file"
  fi
}

export -f extract_taxid

# Get the number of available CPU cores
num_cores=$(nproc)

# Create a temporary file for the filtered database entries
temp_db_file=$(mktemp)

# Extract only the accession.version and taxid whose accession.version starts with NC_
zcat other/nucl_gb.accession2taxid.gz | awk '$2 ~ /^NC_/ {print $2 "\t" $3}' > "$temp_db_file"

# Process each organelle
for organelle in mitochondrion chloroplast; do
  accession_file="${organelle}/other/accessions.txt"
  output_file="${organelle}/other/accessions_taxid.txt"

  echo -e "accession.version\ttaxid" > "$output_file"

  cat "$accession_file" | parallel --gnu -j "$num_cores" extract_taxid {} "$temp_db_file" "$output_file"
done

# Clean up the temporary file
rm "$temp_db_file"