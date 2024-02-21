# Let's start by defining the approach to find the species name for a given taxid
# We will read both the mapping_id_taxa.tsv and taxonomy.tsv files
# Then, for each taxid in mapping_id_taxa.tsv, we will find the corresponding species name in taxonomy.tsv

# File paths for the example
mapping_file_path = 'data/mapping_complete_mt.tsv'
taxonomy_file_path = 'data/taxonomy.tsv'

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

# Now, read the mapping file and find the species name for each taxid
mapping_data = []
with open(mapping_file_path, 'r') as file:
    for line in file:
        fasta_id, taxid = line.strip().split('\t')
        species_name = taxonomy_data.get(taxid, "Unknown")  # Default to "Unknown" if taxid not found
        mapping_data.append((fasta_id, taxid, species_name))

# For demonstration, let's print the results
for item in mapping_data:
    print(item)

# This will output a list of tuples where each tuple contains:
# (fasta identifier, taxid, species name)
# The species name is looked up in the taxonomy_data dictionary using the taxid from the mapping file
