#!/usr/bin/env python3

"""
Description:
This script organizes genome files based on taxonomic criteria provided in the lineage file 
and the ID to species mapping file. It reads the lineage information and the ID to species 
mapping, then organizes the genome files into directories based on genus, family, and order 
taxonomic levels. The directories are created only when 10 or more genomes are present, 
ensuring a threshold for meaningful categorization.

Usage:
./scripts/organize_genomes.py <lineage file> <id and species list> <accession list> <input directory> <output directory>

Arguments:
<lineage file>: Path to the lineage file containing taxonomic lineage information.
<id and species list>: Path to the file containing the mapping of genome IDs to species names.
<accession list>: Path to the file containing the list of accession numbers to sort.
<input directory>: Path to the directory containing input FASTA files.
<output directory>: Path to the base directory where organized files will be stored.

Example usage following PlantEuka folder organization:
python scripts/organize_genomes_after_pw.py other/lineage.tsv other/acc_taxid_name_cp.txt other/accession_list_pw.txt genomes/chloroplast/original genomes/chloroplast/sorted_pw
"""

# Libraries
from collections import defaultdict
import os
import shutil
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 6:
    print("Usage: ./scripts/organize_genomes.py <lineage file> <id and species list> <accession list> <input directory> <output directory>")
    sys.exit(1)

# Extract command-line arguments
lineage_file_path = sys.argv[1]
id_species_file = sys.argv[2]
accession_list_file = sys.argv[3]
input_directory = sys.argv[4]
output_base_dir = sys.argv[5]

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

# Parse accession list
accession_set = set()
with open(accession_list_file, 'r') as f:
    for line in f:
        accession_set.add(line.strip())

# Initialize aggregation structures for organizing FASTA files
species_by_genus = defaultdict(list)
species_by_family = defaultdict(list)
species_by_order = defaultdict(list)

# Aggregate species IDs by their lineage information, only including those in the accession list
for id, species_name in id_to_species.items():
    if species_name in species_to_lineage and id in accession_set:
        genus, family, order = species_to_lineage[species_name]
        species_by_genus[genus].append(id)
        species_by_family[family].append(id)
        species_by_order[order].append(id)

def ensure_directory(dir_path):
    """
    Ensures that the specified directory exists. If the directory does not exist,
    it creates the directory.

    Args:
    dir_path (str): Path to the directory to be ensured.

    Returns:
    None
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def copy_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids):
    """
    Copies FASTA files based on taxonomy criteria if certain conditions are met,
    avoiding duplication across taxonomic levels. It checks if the number of unprocessed
    IDs for a given taxonomy level is sufficient for copying files.

    Args:
    species_dict (defaultdict): Dictionary containing species IDs organized by taxonomy.
    taxonomy (str): Taxonomic level (e.g., genus, family, order).
    processed_ids (set): Set containing processed IDs to prevent duplication.

    Returns:
    None
    """
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
                        shutil.copy(source_file, target_file)  # Use copy instead of move
                        processed_ids.add(id)

# Initialize set of processed IDs to prevent duplication
processed_ids = set()

# Organize and copy files by genus, family, and order, avoiding duplication
for taxonomy, species_dict in [('genus', species_by_genus), ('family', species_by_family), ('order', species_by_order)]:
    copy_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids)
