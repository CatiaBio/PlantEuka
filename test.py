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

# Replace 'your_file.tsv' with the path to your actual file
file_path = 'data/Chlorophyta_genomes/taxon_info.csv'
taxon_id_list = read_taxon_ids(file_path)

# This list will hold your rows
data_to_write = []

for taxon in taxon_id_list:
    try:
        taxon_id = int(taxon.strip())  # Ensure conversion to integer and strip any whitespace
        t = taxoniq.Taxon(taxon_id)  # The taxon id we want to get the genus information
        ranked_lineage = [(t.rank.name, t.scientific_name) for t in t.ranked_lineage]
        genus = str(t.ranked_lineage[1])
        genusTaxID = str(re.findall(r'\d+', genus)[0])
        genus_name = next((name for rank, name in ranked_lineage if rank == 'genus'), None)
        
        # Append each record to the data_to_write list
        data_to_write.append([taxon, genus_name, genusTaxID])

    except KeyError as e:
        genus_name = "NA"
        genusTaxID = "NA"
        data_to_write.append([taxon, genus_name, genusTaxID])
        print(f"TaxonID {taxon} not found in the taxoniq database: {e}")
    except ValueError as e:
        print(f"Invalid TaxonID format: {taxon} - {e}")

# Now write the data to a CSV file
with open('genus_info.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['TaxonID', 'GenusName', 'GenusTaxonID'])  # Writing the header
    writer.writerows(data_to_write)  # Writing the data