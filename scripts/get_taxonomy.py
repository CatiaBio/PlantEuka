import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Process taxonomy files.')
parser.add_argument('--nodes_file', required=True, help='Path to the nodes.dmp file')
parser.add_argument('--names_file', required=True, help='Path to the names.dmp file')
parser.add_argument('--output_file', required=True, help='Path to the output file')

args = parser.parse_args()

# Use the arguments for file paths
nodes_file = args.nodes_file
names_file = args.names_file
output_file = args.output_file

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
