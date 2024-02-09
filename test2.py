import taxoniq
import csv
import re

def read_taxon_ids(file_path):
    taxon_ids = []
    
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')  # Assuming a tab-delimited file
        for row in reader:
            taxon_ids.append(row['TaxonID'])
    
    return taxon_ids

def get_genus_info(taxon_id):
    # Fetches the genus information using the taxoniq library
    t = taxoniq.Taxon(taxon_id)
    for t.rank.name, t.scientific_name in t.ranked_lineage:
        if t.rank.name == 'genus':
            return t.rank.name, t.scientific_name  # Assuming rank.id gives the TaxonID of the genus
    return None, None

# Read the CSV file
file_path = 'data/Chlorophyta_genomes/taxon_info.csv'
rows = read_taxon_ids(file_path)

# Prepare the new CSV data with additional genus information
new_rows = []
for row in rows:
    taxon_id = int(row.strip())  # Convert TaxonID to integer
    genus_name, genus_taxon_id = get_genus_info(taxon_id)  # Get genus info for this taxon_id
    row['Genus'] = genus_name
    row['GenusTaxonID'] = genus_taxon_id
    new_rows.append(row)

# Write the new data to a CSV file
output_file_path = 'data/Chlorophyta_genomes/taxon_info_with_genus.csv'
with open(output_file_path, 'w', newline='') as csvfile:
    fieldnames = list(rows[0].keys())  # Get the fieldnames from the first row
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')
    writer.writeheader()
    writer.writerows(new_rows)

print("Done. The data has been saved to", output_file_path)
