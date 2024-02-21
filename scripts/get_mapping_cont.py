import os
import gzip
import argparse
from Bio import SeqIO

# Setup argument parsing
parser = argparse.ArgumentParser(description='Create a mapping file from FASTA identifiers to NCBI taxonomy identifiers.')
parser.add_argument('--mapping_file', required=True, help='Path to the mapping file containing FASTA identifiers and species names')
parser.add_argument('--taxonomy_file', required=True, help='Path to the taxonomy file containing species names and taxonomy IDs')
parser.add_argument('--output_file', required=True, help='Path to the output file for the new mapping')
args = parser.parse_args()

# Read the mapping file
mapping_data = {}
with open(args.mapping_file, 'r') as file:
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            mapping_data[parts[0]] = parts[1]  # fasta_identifier -> species_name

# Read the taxonomy file
taxonomy_data = {}
with open(args.taxonomy_file, 'r') as file:
    next(file)  # Skip the header
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) > 1:
            species_name = parts[1]  # Name
            taxaid = parts[0]  # Directly use TaxonID
            taxonomy_data[species_name] = taxaid  # species_name -> taxaid

# Replace the species name with the taxaid in the mapping data
new_mapping = [(fasta_id, taxonomy_data.get(species_name, "NA")) for fasta_id, species_name in mapping_data.items()]

# Save the new mapping to the specified output file
with open(args.output_file, 'w') as file:
    for fasta_id, taxaid in new_mapping:
        file.write(f"{fasta_id}\t{taxaid}\n")
