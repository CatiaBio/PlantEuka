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
input_tsv = 'data/viridiplantae_taxonomy_list.tsv'  # Replace with the path to your TSV file
output_txt = 'data/viridiplantae_genus_list.tsv'  # Replace with your desired output file path
extract_genus_names(input_tsv, output_txt)
