#!/usr/bin/env python3

"""
Description:
This script generates a lineage list with additional columns for ID and TaxID based on taxonomic information.
It requires a taxonomy file and a mapping file that associates genome IDs with TaxIDs.

Usage:
./generate_lineage_list.py <taxonomy file> <id_taxid file> <output file>

Arguments:
<taxonomy file>: Path to the taxonomy file (taxonomy.tsv) obtained from the script generate_taxonomy_list.py
<id_taxid file>: Path to the file containing ID to TaxID mappings (id_taxid_cp.tsv)
<output file>: Path for the output file (lineage_with_id.tsv) where the lineage list will be written
"""

import csv
import sys

# Validate command-line arguments
if len(sys.argv) != 4:
    print("Usage: ./generate_lineage_list.py <taxonomy file> <id_taxid file> <output file>")
    sys.exit(1)

taxonomy_file = sys.argv[1]
id_taxid_file = sys.argv[2]
output_file = sys.argv[3]

def load_id_taxid_mapping(id_taxid_file):
    """
    Loads the ID to TaxID mapping from a file.
    """
    mapping = {}
    with open(id_taxid_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                mapping[row[0]] = row[1]
    return mapping

def process_taxonomy_file(taxonomy_file, id_taxid_mapping, output_file):
    """
    Processes the taxonomy file to extract and write taxonomic and ID information.
    """
    lineage_dict = {}
    with open(taxonomy_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                taxon_id, name, rank, lineage = parts
                lineage_dict[taxon_id] = {'rank': rank, 'name': name, 'lineage': lineage.split(',')}

    with open(output_file, 'w', newline='') as f_out:
        csv_writer = csv.writer(f_out, delimiter='\t')
        csv_writer.writerow(['ID', 'TaxID', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum'])

        for taxon_id, info in lineage_dict.items():
            if info['rank'].lower() == 'species':
                rank_info = {'species': info['name'], 'genus': '', 'family': '', 'order': '', 'class': '', 'phylum': ''}
                for ancestor_id in info['lineage']:
                    ancestor_info = lineage_dict.get(ancestor_id, {})
                    ancestor_rank = ancestor_info.get('rank', '').lower()
                    if ancestor_rank in rank_info:
                        rank_info[ancestor_rank] = ancestor_info.get('name', '')
                # Fetch the corresponding ID and TaxID
                for id, taxid in id_taxid_mapping.items():
                    if taxon_id == taxid:
                        csv_writer.writerow([id, taxon_id, rank_info['species'], rank_info['genus'], rank_info['family'], rank_info['order'], rank_info['class'], rank_info['phylum']])

    print(f"Output written to {output_file}")

# Load the ID to TaxID mapping
id_taxid_mapping = load_id_taxid_mapping(id_taxid_file)

# Process the taxonomy file with the loaded mapping
process_taxonomy_file(taxonomy_file, id_taxid_mapping, output_file)