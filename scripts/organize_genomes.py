#!/usr/bin/env python3

"""
Description:
aka sort_genomes.py  

This script organizes genome files based on taxonomic criteria provided in the lineage file 
and the ID to species mapping file. It reads the lineage information and the ID to species 
mapping, then organizes the genome files into directories based on genus, family, and order 
taxonomic levels. The directories are created only when 10 or more genomes are present, 
ensuring a threshold for meaningful categorization. Lists are created for sorted and unsorted 
accession numbers.
"""

# Libraries
from collections import defaultdict
import os
import shutil
import sys

# For Snakemake compatibility
organelle = snakemake.params.organelle

lineage_file_path = snakemake.input.lineage
accession_taxid = snakemake.input.accession_taxid
input_directory = snakemake.input.genomes_dir
output_base_dir = f"{organelle}/genomes/sorted"
sorted_list_path = snakemake.output.sorted_list
unsorted_list_path = snakemake.output.unsorted_list

# Ensure output directories exist
os.makedirs(os.path.dirname(sorted_list_path), exist_ok=True)
os.makedirs(output_base_dir, exist_ok=True)

# Parse full lineage information
species_to_lineage = {}
with open(lineage_file_path, 'r') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            species_to_lineage[parts[0]] = (parts[2], parts[3], parts[4])  # (Genus, Family, Order)

# Parse ID to species mapping
accession_to_taxid = {}
with open(accession_taxid, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            accession_to_taxid[parts[0]] = parts[1]

# Initialize aggregation structures for organizing FASTA files
species_by_genus = defaultdict(list)
species_by_family = defaultdict(list)
species_by_order = defaultdict(list)

# Aggregate species IDs by their lineage information
for accession, taxid in accession_to_taxid.items():
    if taxid in species_to_lineage:
        genus, family, order = species_to_lineage[taxid]
        species_by_genus[genus].append(accession)
        species_by_family[family].append(accession)
        species_by_order[order].append(accession)

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

def copy_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids, sorted_list):
    """
    Copies FASTA files based on taxonomy criteria if certain conditions are met,
    avoiding duplication across taxonomic levels. It checks if the number of unprocessed
    IDs for a given taxonomy level is sufficient for copying files.

    Args:
    species_dict (defaultdict): Dictionary containing species IDs organized by taxonomy.
    taxonomy (str): Taxonomic level (e.g., genus, family, order).
    processed_ids (set): Set containing processed IDs to prevent duplication.
    sorted_list (list): List to store sorted accession numbers.

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
                        shutil.copy(source_file, target_file)
                        processed_ids.add(id)
                        sorted_list.append(id)

# Initialize set of processed IDs to prevent duplication
processed_ids = set()

# Initialize lists for sorted and unsorted accession numbers
sorted_list = []
unsorted_list = []

# Organize and copy files by genus, family, and order, avoiding duplication
for taxonomy, species_dict in [('genus', species_by_genus), ('family', species_by_family), ('order', species_by_order)]:
    copy_fasta_files_if_criteria_met(species_dict, taxonomy, processed_ids, sorted_list)

# Identify unsorted accession numbers
for accession in accession_to_taxid:
    if accession not in sorted_list:
        unsorted_list.append(accession)

# Write sorted and unsorted accession numbers to their respective files
with open(sorted_list_path, 'w') as f:
    for acc in sorted_list:
        f.write(f"{acc}\n")

with open(unsorted_list_path, 'w') as f:
    for acc in unsorted_list:
        f.write(f"{acc}\n")

print("Genome files have been organized and copied successfully.")
print(f"Sorted accession numbers written to {sorted_list_path}")
print(f"Unsorted accession numbers written to {unsorted_list_path}")
