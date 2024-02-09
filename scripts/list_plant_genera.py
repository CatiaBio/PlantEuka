from Bio import Entrez
import csv

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'

def get_plant_genera_taxids_and_names():
    search_term = "txid33090[Subtree] AND (genus[Rank])"
    handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    taxon_ids = record['IdList']
    
    # Open a file to write the results
    with open('plant_genera_taxids_names.tsv', 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        
        # Fetch detailed taxonomic information for each ID
        for taxon_id in taxon_ids:
            handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for record in records:
                genus_name = record.get('ScientificName')
                writer.writerow([genus_name])  # Write each row to the file
    
    print("Data has been saved to plant_genera_taxids_names.tsv")

# Call the function to retrieve and save data
get_plant_genera_taxids_and_names()
