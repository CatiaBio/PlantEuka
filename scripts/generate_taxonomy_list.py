#!/usr/bin/env python3

import sys

# Example usage: ./scripts/generate_taxonomy_list.py other/nodes.dmp other/names.dmp other/taxonomy.tsv
if len(sys.argv) != 4:
    print("Usage: ./generate_taxonomy_list.py <nodes_file> <names_file> <output_file>")
    sys.exit(1)
nodes_file = sys.argv[1]                # path for nodes.dmp file from taxdump.tar.gz
names_file = sys.argv[2]                # path for names.dmp file from taxdump.tar.gz
output_file = sys.argv[3]               # path for output file 

# Read the names and nodes
name_dict = {}
parent_dict = {}
rank_dict = {}

with open(names_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        tax_name = split_line[1].strip()
        name_class = split_line[3].strip()
        if name_class == 'scientific name':
            name_dict[tax_id] = tax_name

with open(nodes_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        parent_tax_id = split_line[1].strip()
        rank = split_line[2].strip()
        parent_dict[tax_id] = parent_tax_id
        rank_dict[tax_id] = rank

# Function to construct the full lineage for a taxon
def get_lineage(tax_id, parent_dict):
    lineage = [tax_id]
    while tax_id in parent_dict and parent_dict[tax_id] != '1':
        tax_id = parent_dict[tax_id]
        lineage.append(tax_id)
    return reversed(lineage)

# Function to check if taxon is in the lineage of Taxon ID 33090
def is_in_lineage_of_33090(tax_id, parent_dict):
    while tax_id in parent_dict and tax_id != '1':
        if tax_id == '33090':
            return True
        tax_id = parent_dict[tax_id]
    return False

# Write the output file
with open(output_file, 'w') as f:
    f.write("TaxonID\tName\tRank\tLineage\n")  # Write the header
    for tax_id in name_dict:
        if is_in_lineage_of_33090(tax_id, parent_dict):
            lineage = get_lineage(tax_id, parent_dict)
            lineage_str = ','.join(lineage)
            rank = rank_dict.get(tax_id, 'no rank')
            f.write(f"{tax_id}\t{name_dict[tax_id]}\t{rank}\t{lineage_str}\n")

print(f"Output written to {output_file}")
