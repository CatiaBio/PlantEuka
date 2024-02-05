from Bio import Entrez
import pandas as pd

# Set your email
Entrez.email = '***REMOVED***'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def read_accession_from_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                data.append(fields[0])
    return data

def fetch_genus_taxid_for_accessions(accession_list):
    data = []
    for accession in accession_list:
        try:
            # Use ESummary to retrieve taxonomy data
            handle = Entrez.esummary(db="nuccore", id=accession)
            record = Entrez.read(handle)
            tax_id = record[0]["TaxId"]
            genus = record[0]["Title"].split()[0]  # Extract genus from title
            
            data.append({"Accession": accession, "Genus": genus, "GenusTaxID": tax_id})
        except Exception as e:
            print(f"Error for accession {accession}: {str(e)}")

    return pd.DataFrame(data)

def save_to_tsv(data_frame, file_name="output.tsv"):
    data_frame.to_csv(file_name, sep='\t', index=False)

# Reading Accession names from your tab-separated file
file_path = "data/chloroplast_taxon_info.tsv"  # Update this path to your actual file location
accession_list = read_accession_from_file(file_path)

# Processing Accession names
genus_taxid_df = fetch_genus_taxid_for_accessions(accession_list)
save_to_tsv(genus_taxid_df, "accessions_to_genus_taxid.tsv")

print("Done. The data has been saved to accessions_to_genus_taxid.tsv")
