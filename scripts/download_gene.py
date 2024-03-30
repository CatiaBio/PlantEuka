#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import sys

# Ensure the correct number of arguments are passed
if len(sys.argv) != 3:
    print("Usage: ./download_genomes.py <query_string> <output_file>")
    sys.exit(1)

query = sys.argv[1]  # Search query string
output_file = sys.argv[2]  # Path for output file

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '92573590fbfb9479ab2167a4f133ee31a408'

def search_genomes(search_term, database='nuccore'):
    """Search GenBank records based on a query."""
    id_list = []
    handle = Entrez.esearch(db=database, term=search_term, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    count = int(record['Count'])
    batch_size = 500
    for start in range(0, count, batch_size):
        handle = Entrez.esearch(db=database, term=search_term, retstart=start, retmax=batch_size)
        record = Entrez.read(handle)
        handle.close()
        id_list.extend(record['IdList'])
    return id_list

def fetch_genome_info(id_list, database='nuccore'):
    """Fetch GenBank records for the given list of IDs."""
    handle = Entrez.efetch(db=database, id=",".join(id_list), rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    return list(records)

def extract_psaA_gene_info(records):
    """Extract the psaA gene start and end positions from records."""
    psaA_gene_info = []
    for record in records:
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                if "psaA" in feature.qualifiers["gene"]:
                    location = str(feature.location).strip("<>")
                    psaA_gene_info.append((record.id, location))
    return psaA_gene_info

def save_psaA_gene_info(psaA_gene_info, output_file_path):
    """Save the psaA gene information to a file."""
    with open(output_file_path, "w") as f:
        for id, location in psaA_gene_info:
            f.write(f"{id} {location}\n")

def main(search_term, output_file_path):
    ids = search_genomes(search_term)
    if ids:
        records = fetch_genome_info(ids)
        psaA_gene_info = extract_psaA_gene_info(records)
        save_psaA_gene_info(psaA_gene_info, output_file_path)
        print(f"psaA gene information saved to {output_file_path}")

if __name__ == "__main__":
    main(query, output_file)
