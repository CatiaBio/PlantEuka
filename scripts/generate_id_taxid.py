#!/usr/bin/env python3

import os
import gzip
import argparse
import sys
from Bio import SeqIO

"""
    Description: 'Create a mapping file from FASTA identifiers to NCBI taxonomy identifiers and include species names.'
    Usage for chloroplast:
    ./scripts/generate_id_taxid.py genomes/chloroplast/id_species_cp.tsv other/taxonomy.tsv genomes/chloroplast/id_taxid_cp.tsv 
    Usage for mitochondrion:
    ./scripts/generate_id_taxid.py genomes/mitochondrion/id_species_mt.tsv other/taxonomy.tsv genomes/mitochondrion/id_taxid_mt.tsv 
"""

# Setup argument parsing
if len(sys.argv) != 4:
    print("Usage: ./scripts/generate_id_taxid.py <id species list> <taxonomy list> <output directory>")
    sys.exit(1)
    
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

