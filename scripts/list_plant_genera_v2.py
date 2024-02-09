from Bio import Entrez
import csv

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'

def search_all_genomes(search_term, database='taxonomy'):
    # Start with an esearch to get the count of total available records
    handle = Entrez.esearch(db=database, term=search_term, retmax=0)
    record = Entrez.read(handle)
    handle.close()

    count = int(record['Count'])  # Total number of records matching the search term
    batch_size = 500  # NCBI recommends not exceeding 500 for large queries
    id_list = []

    # Paginate through the records in batches
    for start in range(0, count, batch_size):
        handle = Entrez.esearch(db=database, term=search_term, retstart=start, retmax=batch_size)
        record = Entrez.read(handle)
        handle.close()
        id_list.extend(record['IdList'])

    return id_list

def get_plant_genera_names(taxon_ids):
    # Open a file to write the results
    with open('data/plant_genera_names.tsv', 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        
        # Fetch detailed taxonomic information for each ID
        for taxon_id in taxon_ids:
            handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for record in records:
                genus_name = record.get('ScientificName')
                writer.writerow([genus_name])  # Write each row to the file
    
    print("Data has been saved to plant_genera_names.tsv")

# Usage
search_term = "txid33090[Subtree] AND (genus[Rank])"
genus_ids = search_all_genomes(search_term, 'taxonomy')
get_plant_genera_names(genus_ids)
