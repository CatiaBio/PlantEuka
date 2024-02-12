import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Process taxonomy files.')
parser.add_argument('--taxonomy_file', required=True, help='Path to the viridiplantae_taxonomy.tsv file')
parser.add_argument('--output_file', required=True, help='Path to the output file')

args = parser.parse_args()

# Use the arguments for file paths
taxonomy_file = args.taxonomy_file
output_file = args.output_file


def extract_genus_names(input_tsv, output_txt):
    genus_names = []

    with open(input_tsv, 'r') as file:
        # Skip header
        next(file)
        
        # Process each line in the TSV file
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:  # Ensure there are enough columns
                taxon_id, name, rank, lineage = parts
                if rank.lower() == 'genus':  # Check if the rank is 'genus'
                    genus_names.append(name)

    # Write the genus names to the output text file
    with open(output_txt, 'w') as file:
        for name in genus_names:
            file.write(name + '\n')

    print(f"Genus names have been saved to {output_txt}")

# Usage
extract_genus_names(taxonomy_file, output_file)
