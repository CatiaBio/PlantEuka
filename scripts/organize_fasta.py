import os
from collections import defaultdict
import shutil

# Corrected file paths assuming the starting point is the current working directory
lineage_file_path = 'data/full_lineage.tsv'
mapping_file_path = 'data/mitochondrion/id_species_mapping.txt'
data_directory = 'data/mitochondrion'  # Base directory for FASTA files

# Parse full lineage
species_to_lineage = {}  # {species_name: (genus, family, order)}
with open(lineage_file_path, 'r') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            species_to_lineage[parts[0]] = (parts[1], parts[2], parts[3])  # (Genus, Family, Order)

# Parse ID to species mapping
id_to_species = {}  # {id: species_name}
with open(mapping_file_path, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            id_to_species[parts[0]] = parts[1]

# Initialize aggregation structures
species_by_genus = defaultdict(list)
species_by_family = defaultdict(list)
species_by_order = defaultdict(list)

# Processed IDs to prevent duplication
processed_ids = set()

# Aggregate species IDs by genus, family, and order based on the species name from id_species_mapping
for id, species_name in id_to_species.items():
    if species_name in species_to_lineage:
        genus, family, order = species_to_lineage[species_name]
        species_by_genus[genus].append(id)
        species_by_family[family].append(id)
        species_by_order[order].append(id)

def ensure_directory(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

# Function to move/copy FASTA files without duplication
def move_fasta_files(id_list, target_dir, processed_ids):
    for id in id_list:
        if id not in processed_ids:
            source_file = os.path.join(data_directory, f"{id}.fasta.gz")
            if os.path.exists(source_file):
                target_file = os.path.join(target_dir, f"{id}.fasta.gz")
                shutil.copy(source_file, target_file)  # Use shutil.move to move instead
                processed_ids.add(id)

# Directory for saving files based on counts
output_base_dir = 'data/mitochondrion'

# Function to count valid FASTA files for each taxonomy
def count_valid_files(id_list):
    count = 0
    for id in id_list:
        if os.path.exists(os.path.join(data_directory, f"{id}.fasta.gz")):
            count += 1
    return count

# Function to move/copy FASTA files without duplication, adjusted to check counts before proceeding
def move_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids):
    for taxon, id_list in species_dict.items():
        # Check if this taxon already has processed IDs to avoid duplication
        unprocessed_ids = [id for id in id_list if id not in processed_ids]
        if unprocessed_ids:
            valid_file_count = count_valid_files(unprocessed_ids)
            if valid_file_count >= 10:
                target_dir = os.path.join(output_base_dir, taxonomy, taxon)
                ensure_directory(target_dir)
                for id in unprocessed_ids:
                    source_file = os.path.join(data_directory, f"{id}.fasta.gz")
                    if os.path.exists(source_file):
                        target_file = os.path.join(target_dir, f"{id}.fasta.gz")
                        shutil.copy(source_file, target_file)  # Use shutil.move to move instead
                        processed_ids.add(id)

# Initialize processed IDs to ensure no duplication across taxonomic levels
processed_ids = set()

# Process and organize files by genus, family, and order without duplication
for taxonomy, species_dict in [('genus', species_by_genus), ('family', species_by_family), ('order', species_by_order)]:
    move_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids)