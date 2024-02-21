import os
import gzip
import argparse
from Bio import SeqIO

# Setup argument parsing
parser = argparse.ArgumentParser(description='Create a mapping file from FASTA identifiers to NCBI taxonomy identifiers and include species names.')
parser.add_argument('--mapping_file', required=True, help='Path to the mapping file containing FASTA identifiers and species names')
parser.add_argument('--taxonomy_file', required=True, help='Path to the taxonomy file containing species names and taxonomy IDs')
parser.add_argument('--output_file', required=True, help='Path to the output file for the new mapping')
args = parser.parse_args()

# Read the mapping file to create a dictionary mapping FASTA identifiers to species names
mapping_data = {}
with open(args.mapping_file, 'r') as file:
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            mapping_data[parts[0]] = parts[1]  # fasta_identifier -> species_name

# Read the taxonomy file to create a dictionary mapping species names to taxonomy IDs
taxonomy_data = {}
with open(args.taxonomy_file, 'r') as file:
    next(file)  # Skip the header
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) > 1:
            species_name = parts[1]  # Name
            taxaid = parts[0]  # TaxonID
            taxonomy_data[species_name] = taxaid  # species_name -> taxaid

# Construct new mapping with FASTA ID, species name, and corresponding taxaid
new_mapping = [(fasta_id, species_name, taxonomy_data.get(species_name, "NA")) for fasta_id, species_name in mapping_data.items()]

# Save the new mapping to the output file, including species names
with open(args.output_file, 'w') as file:
    for fasta_id, species_name, taxaid in new_mapping:
        file.write(f"{fasta_id}\t{species_name}\t{taxaid}\n")
