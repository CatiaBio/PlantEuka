# File paths for the example
mapping_file_path = 'data/mapping_complete_mt.tsv'
taxonomy_file_path = 'data/taxonomy.tsv'
output_file_path = 'data/fasta_species_mapping.tsv'  # Path for the new output TSV file

# Read the taxonomy file to create a dictionary mapping TaxonID to Name
taxonomy_data = {}
with open(taxonomy_file_path, 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        parts = line.strip().split('\t')
        if len(parts) > 1:
            taxon_id = parts[0]  # TaxonID
            name = parts[1]  # Name
            taxonomy_data[taxon_id] = name

# Now, read the mapping file and find the species name for each taxid, then save the result
with open(mapping_file_path, 'r') as file, open(output_file_path, 'w') as output_file:
    for line in file:
        fasta_id, taxid = line.strip().split('\t')
        species_name = taxonomy_data.get(taxid, "Unknown")  # Default to "Unknown" if taxid not found
        output_file.write(f"{fasta_id}\t{species_name}\n")  # Write fasta_id and species_name to output file
