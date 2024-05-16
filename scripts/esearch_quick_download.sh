#!/bin/bash

# File containing the list of parts or IDs, assumed one per line
input_file= "other/part_aa"

# Loop through each line in the input file
while IFS= read -r id
do
    # Use the ID to fetch the sequence and save it to a file named after the ID
    esearch -db nucleotide -query "$id" | efetch -format fasta > "genomes/2405_chloroplast/${id}.fasta"
done < "$input_file"
