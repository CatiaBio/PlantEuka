import os
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
with open('viridiplantae_taxonomy.tsv', 'r') as file:
    next(file)  # Skip the header
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) >= 6:  # Ensure there are enough columns
            species, genus, family, order, class_, phylum = parts
            genus_to_family[genus] = family

# Create directories for families and move genus directories into them
for genus, family in genus_to_family.items():
    family_dir = os.path.join(family_base_dir, family)
    genus_dir = os.path.join(genus_base_dir, genus)
    target_dir = os.path.join(family_dir, genus)

    # Create the family directory if it doesn't exist
    os.makedirs(family_dir, exist_ok=True)

    # Move the genus directory into the family directory (uncomment the desired operation)
    # shutil.move(genus_dir, target_dir)                        # Use this to move directories
    shutil.copytree(genus_dir, target_dir, dirs_exist_ok=True)  # Use this to copy directories
