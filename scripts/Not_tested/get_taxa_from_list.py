from Bio import Entrez
import pandas as pd

# Set your email
Entrez.email = '***REMOVED***'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def read_ids_from_file(file_path):
    with open(file_path, 'r') as file:
        ids = [line.strip() for line in file if line.strip()]
    return ids

def fetch_genus_taxid_for_ids(id_list):
    data = []
    for id in id_list:
        # NCBI ELink to find related taxonomy IDs
        handle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=id, linkname="nuccore_taxonomy")
        link_results = Entrez.read(handle)
        if link_results[0]["LinkSetDb"]:
            tax_ids = [link['Id'] for link in link_results[0]['LinkSetDb'][0]['Link']]
            if tax_ids:
                # Fetch taxonomy records
                tax_handle = Entrez.efetch(db="taxonomy", id=tax_ids, retmode="xml")
                tax_records = Entrez.read(tax_handle)
                # Find the genus level information
                for tax_record in tax_records:
                    genus_info = next((line for line in tax_record['LineageEx'] if line['Rank'] == 'genus'), None)
                    if genus_info:
                        data.append({"ID": id, "GenusTaxID": genus_info['TaxId']})
                        break
    return pd.DataFrame(data)

def save_to_tsv(data_frame, file_name="output.tsv"):
    data_frame.to_csv(file_name, sep='\t', index=False)

# Reading IDs from a file
file_path = "id_list.txt"  # Update this path to your actual file location
id_list = read_ids_from_file(file_path)

# Processing IDs
genus_taxid_df = fetch_genus_taxid_for_ids(id_list)
save_to_tsv(genus_taxid_df, "ids_to_genus_taxid.tsv")

print("Done. The data has been saved to ids_to_genus_taxid.tsv")


