#!/usr/bin/env python3

from collections import defaultdict
import os
import shutil
import sys

"""
    Usage for chloroplast:
    ./scripts/organize_genomes.py other/lineage.tsv genomes/chloroplast/id_species_cp.tsv genomes/chloroplast/unsorted genomes/chloroplast/sorted
    Usage for mitochondrion:
    ./scripts/organize_genomes.py other/lineage.tsv genomes/mitochondrion/id_species_mt.tsv genomes/mitochondrion/unsorted genomes/mitochondrion/sorted
"""
# Setup argument parsing
if len(sys.argv) != 5:
    print("Usage: ./scripts/organize_genomes.py <lineage file> <id species list> <input directory> <output directory>")
    sys.exit(1)

lineage_file_path = sys.argv[1] # Path to lineage file 
id_species_file = sys.argv[2]   # Path do id_species file  
input_directory = sys.argv[3]   # Input directory for FASTA files
output_base_dir = sys.argv[4]   # Base directory for organized files

# Parse full lineage information
species_to_lineage = {}
with open(lineage_file_path, 'r') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            species_to_lineage[parts[0]] = (parts[1], parts[2], parts[3])  # (Genus, Family, Order)

# Parse ID to species mapping
id_to_species = {}
with open(id_species_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            id_to_species[parts[0]] = parts[1]

# Initialize aggregation structures for organizing FASTA files
species_by_genus = defaultdict(list)
species_by_family = defaultdict(list)
species_by_order = defaultdict(list)

# Aggregate species IDs by their lineage information
for id, species_name in id_to_species.items():
    if species_name in species_to_lineage:
        genus, family, order = species_to_lineage[species_name]
        species_by_genus[genus].append(id)
        species_by_family[family].append(id)
        species_by_order[order].append(id)

# Ensure the directory exists before trying to move files into it
def ensure_directory(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

# Function to move FASTA files based on taxonomy criteria, ensuring no duplication across levels
def move_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids):
    for taxon, id_list in species_dict.items():
        unprocessed_ids = [id for id in id_list if id not in processed_ids]
        if unprocessed_ids:
            valid_file_count = len(unprocessed_ids)
            if valid_file_count >= 10:
                target_dir = os.path.join(output_base_dir, taxonomy, taxon)
                ensure_directory(target_dir)
                for id in unprocessed_ids:
                    source_file = os.path.join(input_directory, f"{id}.fasta.gz")
                    if os.path.exists(source_file):
                        target_file = os.path.join(target_dir, f"{id}.fasta.gz")
                        shutil.move(source_file, target_file)  # Corrected to move
                        processed_ids.add(id)

# Initialize set of processed IDs to prevent duplication
processed_ids = set()

# Organize and move files by genus, family, and order, avoiding duplication
for taxonomy, species_dict in [('genus', species_by_genus), ('family', species_by_family), ('order', species_by_order)]:
    move_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids)
