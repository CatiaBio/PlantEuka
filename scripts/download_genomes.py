#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import os
import sys

# Example usage: 
#./scripts/download_genomes.py genomes/chloroplast/unsigned/id_species_cp.tsv genomes/chloroplast/unsigned
#./scripts/download_genomes.py genomes/mitochondrion/unsigned/id_species_mt.tsv genomes/mitochondrion/unsigned

if len(sys.argv) != 4:
    print("Usage: ./download_genomes.py <query_string> <id_specie_list> <output_file>")
    sys.exit(1)
query = sys.argv[1]                     # search query string
id_species = sys.argv[2]           # path for id_species file 
output_file = sys.argv[3]               # path for output file 

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def search_genomes(search_term, database='nuccore'):
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
    # Fetch GenBank records for the given list of IDs
    handle = Entrez.efetch(db=database, id=",".join(id_list), rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    return list(records)

def save_id_species_mapping(id_species_mapping, mapping_file_path):
    with open(mapping_file_path, "w") as f:
        for id, species in id_species_mapping:
            f.write(f"{id}\t{species}\n")

def save_genome_sequences(records, directory, mapping_file_path):
    id_species_mapping = []
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    for record in records:
        # fasta_filepath = os.path.join(directory, f"{record.id}.fasta.gz")
        # with gzip.open(fasta_filepath, "wt") as output_handle:
        #     SeqIO.write(record, output_handle, "fasta")
        
        species_name = record.annotations.get('organism', 'Unknown')
        id_species_mapping.append((record.id, species_name))
    
    save_id_species_mapping(id_species_mapping, mapping_file_path)

def main(search_term, output_dir, mapping_file_path):
    ids = search_genomes(search_term)
    if ids:
        records = fetch_genome_info(ids)
        save_genome_sequences(records, output_dir, mapping_file_path)

if __name__ == "__main__":
    main(query, output_file, id_species)
