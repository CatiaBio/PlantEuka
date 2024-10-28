#!/usr/bin/env python3

"""
Description:
The script is designed to generate a taxonomy list based on information from two files 
(nodes_file and names_file) obtained from the NCBI Taxonomy database. 
It extracts taxonomic information, including the taxon ID, name, rank, and lineage, and 
writes it to an output file, taxonomy.tsv. It also processes this taxonomy file 
to generate a lineage list based on taxonomic information, including species, genus, 
family, order, class, and phylum, and writes it to another output file, lineage.tsv.
"""

# Libraries 
import sys
import subprocess
import os
import csv

# Extract command-line arguments
nodes_file = "other/nodes.dmp"  # file to be extracted from taxdump.tar.gz
names_file = "other/names.dmp"  # file to be extracted from taxdump.tar.gz
selected_tax_id = "33090"       # tax ID corresponding to green plants
taxonomy_output_file = "other/taxonomy.tsv"
lineage_output_file = "other/lineage.tsv"

# Path to the tar.gz file
#tar_file = "other/taxdump.tar.gz"

# # Extract the required files from the tar.gz archive
# print("Extracting nodes.dmp and names.dmp from taxdump.tar.gz...")
# subprocess.run(["tar", "-zxvf", tar_file, "-C", "other", "nodes.dmp", "names.dmp"], check=True)

# # Check if the files exist after extraction
# if not os.path.isfile(nodes_file) or not os.path.isfile(names_file):
#     print("Error: nodes.dmp or names.dmp not found after extraction.")
#     sys.exit(1)

# Create dictionaries from names and nodes file
name_dict = {}  # stores scientific names of taxa
parent_dict = {}  # stores parent-child relationships
rank_dict = {}  # stores ranks of taxa 

# Parsing through names_file 
print("Parsing names.dmp...")
with open(names_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        tax_name = split_line[1].strip()
        name_class = split_line[3].strip()
        if name_class == 'scientific name':
            name_dict[tax_id] = tax_name

# Parsing through nodes_file
print("Parsing nodes.dmp...")
with open(nodes_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        parent_tax_id = split_line[1].strip()
        rank = split_line[2].strip()
        parent_dict[tax_id] = parent_tax_id
        rank_dict[tax_id] = rank

# Construct the full lineage for a taxon
def get_lineage(tax_id, parent_dict):
    """
    This function aims to retrieve the complete lineage of a specified taxon 
    ID by traversing through the parent-child relationships provided in the 
    parent dictionary.

    Args:
    tax_id (str): Taxon ID for which lineage is to be retrieved.
    parent_dict (dict): Dictionary mapping taxon IDs to their parent taxon IDs.

    Returns:
    list: List containing the complete lineage of the specified taxon ID.
    """
    lineage = [tax_id]
    while tax_id in parent_dict and parent_dict[tax_id] != '1':
        tax_id = parent_dict[tax_id]
        lineage.append(tax_id)
    return list(reversed(lineage))

# Check if taxon is in the lineage of a selected Taxon ID 
def is_in_lineage_of_x(tax_id, parent_dict, selected_tax_id):
    """
    This function aims to determine whether a specified taxon ID is within the lineage 
    of specific Taxon ID by traversing through the parent-child relationships provided 
    in the parent dictionary.

    Args:
    tax_id (str): Taxon ID to be checked for lineage inclusion.
    parent_dict (dict): Dictionary mapping taxon IDs to their parent taxon IDs.

    Returns:
    bool: True if the specified taxon ID is within the lineage of selected taxon ID,
    False otherwise.
    """
    while tax_id in parent_dict and tax_id != '1':
        if tax_id == selected_tax_id:
            return True
        tax_id = parent_dict[tax_id]
    return False

# Write the taxonomy output file
print("Writing taxonomy.tsv...")
with open(taxonomy_output_file, 'w') as f:
    # Write header with TaxonID, Name, Rank and Lineage 
    f.write("taxid\tname\trank\tlineage\n")    
    for tax_id in name_dict:
        if is_in_lineage_of_x(tax_id, parent_dict, selected_tax_id):
            lineage = get_lineage(tax_id, parent_dict)
            lineage_str = ','.join(lineage)
            rank = rank_dict.get(tax_id, 'no rank')
            f.write(f"{tax_id}\t{name_dict[tax_id]}\t{rank}\t{lineage_str}\n")

print(f"Taxonomy output written to {taxonomy_output_file}")

# Remove the extracted files
os.remove(nodes_file)
os.remove(names_file)

# Function to process the taxonomy file and extract rank information
def process_taxonomy_file(taxonomy_file, output_file):
    """
    This function processes the taxonomy file to extract taxonomic rank information
    for each species. It reads the taxonomy file line by line, extracts the taxon ID,
    name, rank, and lineage, and writes the extracted rank information to the output file.

    Args:
    taxonomy_file (str): Path to the taxonomy file containing taxonomic information.
    output_file (str): Path for the output file where the rank information will be written.

    Returns:
    None
    """
    # Initialize a dictionary to hold taxon ID to its lineage mapping
    lineage_dict = {}
    
    # Open the taxonomy file and read line by line
    with open(taxonomy_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # Ensure the line has enough parts to process
            if len(parts) >= 4:
                taxon_id, name, rank, lineage = parts
                # Store the lineage information, which is a comma-separated list of ancestor taxon IDs
                lineage_dict[taxon_id] = {'rank': rank, 'name': name, 'lineage': lineage.split(',')}

    # Open the output file for writing
    with open(output_file, 'w', newline='') as f_out:
        csv_writer = csv.writer(f_out, delimiter='\t')
        # Write the header row
        csv_writer.writerow(['taxid', 'species', 'genus', 'family', 'order', 'class', 'phylum'])
        
        # Process each taxon ID to find and extract the specified ranks
        for taxon_id, info in lineage_dict.items():
            if info['rank'].lower() == 'species':
                # Initialize a dictionary to hold the rank information for the current taxon
                rank_info = {'taxid': taxon_id, 'species': info['name'], 'genus': '', 'family': '', 'order': '', 'class': '', 'phylum': ''}
                # Traverse the lineage to fill in the rank information
                for ancestor_id in info['lineage']:
                    ancestor_info = lineage_dict.get(ancestor_id, {})
                    ancestor_rank = ancestor_info.get('rank', '').lower()
                    if ancestor_rank in rank_info:
                        rank_info[ancestor_rank] = ancestor_info.get('name', '')
                # Write the extracted rank information for the current species
                csv_writer.writerow([rank_info['taxid'], rank_info['species'], rank_info['genus'], rank_info['family'], rank_info['order'], rank_info['class'], rank_info['phylum']])

    print(f"Lineage output written to {output_file}")

# Run the function with the provided arguments
process_taxonomy_file(taxonomy_output_file, lineage_output_file)