#!/bin/bash

# Define the input and output files
input_file="other/2405_all_mt_genomes.txt"
db_file="other/NC_accession_taxid.txt"
output_file="other/2405_id_taxid_mt.txt"

# Write header to the output file
echo -e "accession.version\ttaxid" > "$output_file"

# Read each accession.version from the input file
while read -r accession
do
    if [[ -n "$accession" ]]; then
        # Use grep to search the database file and extract the line
        line=$(grep "^$accession" "$db_file")

        if [[ -n "$line" ]]; then
            # Extract the taxid from the line, assuming tab as the delimiter
            taxid=$(echo "$line" | cut -f2)
        else
            # If no line is found, set taxid as 'Not_Found'
            taxid="Not_Found"
        fi

        # Write the result to the output file
        echo -e "$accession\t$taxid" >> "$output_file"
    fi
done < "$input_file"