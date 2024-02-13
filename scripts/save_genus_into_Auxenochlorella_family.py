import os
import shutil

# Path to the TSV file containing the taxonomy information
taxonomy_file = 'data/full_lineage.tsv'

# The base directory where the current genus directories are located
genus_base_dir = 'data/chloroplast'

# The target base directory where family directories will be created
family_base_dir = 'data/chloroplast/family'

# Read the TSV file to create a mapping of genus to family
genus_to_family = {}
with open(taxonomy_file, 'r') as file:
    next(file)  # Skip the header
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) >= 6:  # Ensure there are enough columns
            species, genus, family, order, class_, phylum = parts
            genus_to_family[genus] = family

# Filter to process only the genera belonging to a specific family (e.g., Auxenochlorella)
specific_family = "Auxenochlorella"
specific_genus_to_process = {genus: family for genus, family in genus_to_family.items() if family == specific_family}

# Create directories for the specific family and copy genus directories into them
for genus, family in specific_genus_to_process.items():
    family_dir = os.path.join(family_base_dir, family)
    genus_dir = os.path.join(genus_base_dir, genus)
    target_dir = os.path.join(family_dir, genus)

    # Create the family directory if it doesn't exist
    os.makedirs(family_dir, exist_ok=True)

    # Copy the genus directory into the family directory
    shutil.copytree(genus_dir, target_dir, dirs_exist_ok=True)  # Use this to copy directories

print(f"Processed directories for the family '{specific_family}'.")
