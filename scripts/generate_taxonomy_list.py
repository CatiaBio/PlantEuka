#!/usr/bin/env python3

"""
Description:
The script is designed to generate a taxonomy list based on information from two files 
(nodes_file and names_file) obtained from the NCBI Taxonomy database. 
It extracts taxonomic information, including the taxon ID, name, rank, and lineage, and 
writes it to an output file specified by the user.

Usage: 
./generate_taxonomy_list.py <taxdump nodes file> <taxdump names file> <output file>

Arguments:
<taxdump nodes file>: Path to the nodes.dmp file obtained from the NCBI Taxonomy database.
<taxdump names file>: Path to the names.dmp file obtained from the NCBI Taxonomy database.
<output file>: Path for the output file where the taxonomy list will be written.

Example usage following PlantEuka folder organization:
./scripts/generate_taxonomy_list.py other/taxdump/nodes.dmp other/taxdump/names.dmp other/taxonomy.tsv
"""

# Libraries 
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 4:
    print("Usage: ./generate_taxonomy_list.py <taxdump nodes file> <taxdump names file> <output file>")
    sys.exit(1)

# Extract command-line arguments
nodes_file = sys.argv[1]                
names_file = sys.argv[2]                
output_file = sys.argv[3]               

# Create dictionaries from names and nodes file
name_dict = {}                                   # stores scientific names of taxa
parent_dict = {}                                 # stores parent-child relationships
rank_dict = {}                                   # stores ranks of taxa 

# Parsing trough names_file 
with open(names_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        tax_name = split_line[1].strip()
        name_class = split_line[3].strip()
        if name_class == 'scientific name':
            name_dict[tax_id] = tax_name

# Parsing trough nodes_file
with open(nodes_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        parent_tax_id = split_line[1].strip()
        rank = split_line[2].strip()
        parent_dict[tax_id] = parent_tax_id
        rank_dict[tax_id] = rank

# Cnstruct the full lineage for a taxon
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
    return reversed(lineage)

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

# Write the output file
with open(output_file, 'w') as f:
    # Write header with TaxonID, Name, Rank and Lineage 
    f.write("TaxonID\tName\tRank\tLineage\n")    
    for tax_id in name_dict:
        # Specify the tax id to create the list for 
        selected_tax_id = 33090 # Taxon ID 33090 represents green plants
        if is_in_lineage_of_x(tax_id, parent_dict, selected_tax_id):
            lineage = get_lineage(tax_id, parent_dict,selected_tax_id)
            lineage_str = ','.join(lineage)
            rank = rank_dict.get(tax_id, 'no rank')
            f.write(f"{tax_id}\t{name_dict[tax_id]}\t{rank}\t{lineage_str}\n")

print(f"Output written to {output_file}")