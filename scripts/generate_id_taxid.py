#!/usr/bin/env python3

"""
Description: 
This script generates a mapping file from FASTA identifiers to NCBI taxonomy identifiers 
and includes species names. It reads the ID to species mapping file and the taxonomy 
list, then creates a new mapping file containing FASTA identifiers, species names, 
and corresponding taxonomy identifiers.

Usage:
./generate_id_taxid.py <id species list> <taxonomy list> <output directory>

Arguments:
<id species list>: Path to the file containing the mapping of FASTA identifiers to species names.
<taxonomy list>: Path to the taxonomy list file containing species names and taxonomy IDs.
<output directory>: Path for the output file where the new mapping will be saved.

Example usage following PlantEuka folder organization:
./scripts/generate_id_taxid.py genomes/chloroplast/id_species_cp.tsv other/taxonomy.tsv genomes/chloroplast/id_taxid_cp.tsv 
"""

# Libraries 
import sys
from Bio import SeqIO

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 4:
    print("Usage: ./scripts/generate_id_taxid.py <id species list> <taxonomy list> <output directory>")
    sys.exit(1)

# Extract command-line arguments
mapping_file = sys.argv[1]
taxonomy_file = sys.argv[2]
output_file = sys.argv[3]

# Read the mapping file to create a dictionary mapping FASTA identifiers to species names
mapping_data = {}
with open(mapping_file, 'r') as file:
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            mapping_data[parts[0]] = parts[1]  # fasta_identifier -> species_name

# Read the taxonomy file to create a dictionary mapping species names to taxonomy IDs
taxonomy_data = {}
with open(taxonomy_file, 'r') as file:
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
with open(output_file, 'w') as file:
    for fasta_id, species_name, taxaid in new_mapping:
        file.write(f"{fasta_id}\t{taxaid}\n")