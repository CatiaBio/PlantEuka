# Libraries
from Bio import Entrez
import os

# Access the file paths from Snakemake's input and output
input_file = snakemake.input[0]  
output_file = snakemake.output[0]    

# Read email and API key from the input file
with open(input_file, "r") as info:
    get_info = info.readlines()
    email = get_info[0].strip()
    api_key = get_info[1].strip()

# Set email and API key for Entrez
Entrez.email = email
Entrez.api_key = api_key

# Define the search query
query = "plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]"

def fetch_all_accessions(batch_size=1000):
    accession_numbers = set()
    start = 0
    while True:
        # Perform a paginated search
        with Entrez.esearch(db="nucleotide", term=query, retstart=start, retmax=batch_size, idtype="acc") as search_handle:
            search_results = Entrez.read(search_handle)
            batch_accessions = search_results["IdList"]
            accession_numbers.update(batch_accessions)
            
            # If fewer results than the batch size were returned, we've fetched everything
            if len(batch_accessions) < batch_size:
                break
            
            # Otherwise, move to the next batch
            start += batch_size
    
    return accession_numbers

def read_existing_accessions(file_path):
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            return set(line.strip() for line in file)
    else:
        return set()  # Return an empty set if the file doesn't exist

def save_accessions(file_path, accessions):
    with open(file_path, "w") as file:
        for accession in sorted(accessions):
            file.write(accession + "\n")

# Fetch all accessions in batches and read existing ones if the file exists
latest_accessions = fetch_all_accessions()
existing_accessions = read_existing_accessions(output_file)

# Combine the two sets to get an updated list
updated_accessions = existing_accessions | latest_accessions

# Save the updated list if there are new additions
if updated_accessions != existing_accessions:
    save_accessions(output_file, updated_accessions)
    print(f"Accession list updated with {len(updated_accessions) - len(existing_accessions)} new accession(s).")
else:
    print("No new accessions found.")