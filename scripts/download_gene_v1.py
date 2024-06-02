#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import sys

# Ensure the correct number of arguments are passed
if len(sys.argv) != 3:
    print("Usage: ./download_genomes.py <query_string> <output_file>")
    sys.exit(1)

query = sys.argv[1]  # Search query string
output_file = sys.argv[2]  # Path for output file

# Set your NCBI Entrez email and API key here
Entrez.email = '***REMOVED***'  # Use your actual email
Entrez.api_key = '***REMOVED***'  # Use your actual API key


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

def extract_rbcL_gene_info(records):
    """Extract the rbcL gene start and end positions from records."""
    rbcL_gene_info = []
    for record in records:
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                if "rbcL" in feature.qualifiers["gene"]:
                    location = str(feature.location).strip("<>")
                    rbcL_gene_info.append((record.id, location))
    return rbcL_gene_info

def save_rbcL_gene_info(rbcL_gene_info, output_file_path):
    """Save the rbcL gene information to a file."""
    with open(output_file_path, "a") as f:
        for id, location in rbcL_gene_info:
            f.write(f"{id}\t{location}\n")

def main(search_term, output_file_path):
    ids = search_genomes(search_term)
    if ids:
        for i in range(0, len(ids), 100):  # Fetch in chunks to avoid overloading NCBI
            batch_ids = ids[i:i + 100]
            records = fetch_genome_info(batch_ids)
            rbcL_gene_info = extract_rbcL_gene_info(records)
            save_rbcL_gene_info(rbcL_gene_info, output_file_path)
            print(f"Processed {len(batch_ids)} IDs, total processed: {i + len(batch_ids)}")

if __name__ == "__main__":
    main(query, output_file)
