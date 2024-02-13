import argparse
import csv

# Setup argument parsing
parser = argparse.ArgumentParser(description='Extract and organize taxonomy ranks from lineage.')
parser.add_argument('--taxonomy_file', required=True, help='Path to the taxonomy file')
parser.add_argument('--output_file', required=True, help='Path to save the output file')
args = parser.parse_args()

# Function to process the taxonomy file and extract rank information
def process_taxonomy_file(taxonomy_file, output_file):
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
        csv_writer.writerow(['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum'])
        
        # Process each taxon ID to find and extract the specified ranks
        for taxon_id, info in lineage_dict.items():
            if info['rank'].lower() == 'species':
                # Initialize a dictionary to hold the rank information for the current taxon
                rank_info = {'species': info['name'], 'genus': '', 'family': '', 'order': '', 'class': '', 'phylum': ''}
                # Traverse the lineage to fill in the rank information
                for ancestor_id in info['lineage']:
                    ancestor_info = lineage_dict.get(ancestor_id, {})
                    ancestor_rank = ancestor_info.get('rank', '').lower()
                    if ancestor_rank in rank_info:
                        rank_info[ancestor_rank] = ancestor_info.get('name', '')
                # Write the extracted rank information for the current species
                csv_writer.writerow([rank_info['species'], rank_info['genus'], rank_info['family'], rank_info['order'], rank_info['class'], rank_info['phylum']])

    print(f"Output written to {output_file}")

# Run the function with the provided arguments
process_taxonomy_file(args.taxonomy_file, args.output_file)
