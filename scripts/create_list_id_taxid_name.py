# Define the file paths
taxonomy_file = 'other/2405_taxid.tsv'
accession_file = 'other/2405_id_taxid_mt.txt'
output_file = 'other/2405_id_taxid_name_mt.txt'

# Read taxonomy information into a dictionary
taxid_to_name = {}
with open(taxonomy_file, 'r') as f:
    # Skip the header
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) > 2:
            taxid = parts[0]
            name = parts[1]
            taxid_to_name[taxid] = name

# Read accession numbers and their corresponding TaxonIDs and write the results
with open(accession_file, 'r') as f, open(output_file, 'w') as out:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            accession = parts[0]
            taxid = parts[1]
            # Lookup the name using the TaxonID
            name = taxid_to_name.get(taxid, 'Name_Not_Found')
            # Write to the output file
            out.write(f'{accession}\t{taxid}\t{name}\n')

print("Output has been written to", output_file)
