#!/bin/bash

# Define the input and output files
input_file="other/2405_all_cp_genomes.txt"
db_file="NC_accession_taxid.txt"
output_file="other/2405_id_taxid_cp.txt"

# Write header to the output file
echo -e "accession.version\ttaxid" > "$output_file"

# Read each accession.version from the input file
while read -r accession
do
    if [[ -n "$accession" ]]; then
        # Use zgrep to search the gzipped database file and extract the line
        line=$(zgrep "^$accession\t" "$db_file")

        if [[ -n "$line" ]]; then
            # Extract the taxid from the line
            taxid=$(echo "$line" | cut -f3)
        else
            # If no line is found, set taxid as 'Not_Found'
            taxid="Not_Found"
        fi

        # Write the result to the output file
        echo -e "$accession\t$taxid" >> "$output_file"
    fi
done < "$input_file"
