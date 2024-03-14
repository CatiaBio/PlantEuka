import os
from collections import defaultdict
import shutil
import argparse

"""
    Usage for chloroplast:
    python organize_fasta.py --lineage_file ../data/full_lineage.tsv --mapping_file ../data/mapping_id_species_cp.txt --input_dir ../data/chloroplast/raw --output_dir ../data/chloroplast/organized
    Usage for mitochondrion:
    python organize_fasta.py --lineage_file ../data/full_lineage.tsv --mapping_file ../data/mapping_id_species_mt.txt --input_dir ../data/mitochondrion/raw --output_dir ../data/mitochondrion/organized
"""
# Setup argument parsing
parser = argparse.ArgumentParser(description='Organize sequence files per genus, family or order.')
parser.add_argument('--lineage_file', required=True, help='Path to the lineage file')
parser.add_argument('--mapping_file', required=True, help='Path to the mapping file')
parser.add_argument('--input_dir', required=True, help='Path to the input directory containing FASTA files')
parser.add_argument('--output_dir', required=True, help='Path to the output directory')
args = parser.parse_args()

# Adjusted file paths based on argument parsing
lineage_file_path = args.lineage_file
mapping_file_path = args.mapping_file
input_directory = args.input_dir  # Changed to reflect the actual input directory for FASTA files
output_base_dir = args.output_dir  # This will be used as the base directory for organized files

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
with open(mapping_file_path, 'r') as f:
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
