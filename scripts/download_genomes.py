#!/usr/bin/env python3

"""
Description:
This script is designed to download genome sequences from the NCBI database based on a provided search query. 
It retrieves genome records matching the specified search query, saves the sequences in FASTA format, 
and generates a mapping file linking genome IDs to species names.

Usage:
./download_genomes.py <query string> <id and species list> <output directory>

Arguments:
<query string>: NCBI search query string used to retrieve genome sequences matching specific criteria.
<id and species list>: Path to the file where the mapping of genome IDs to species names will be saved (e.g., id_species_cp.tsv).
<output folder>: Path to the directory where the downloaded genome sequences will be stored in FASTA format.

Example: 
./scripts/download_genomes.py "plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]" genomes/chloroplast/id_species_cp.tsv genomes/chloroplast/unsorted
"""

# Libraries
from Bio import Entrez, SeqIO
import os
import sys
import gzip

# Check if the correct number of command-line arguments are provided. 
if len(sys.argv) != 4:
    print("Usage: ./download_genomes.py <query string> <id and species list> <output directory>")
    sys.exit(1)

# Extract command-line arguments
query = sys.argv[1]                   
id_species = sys.argv[2]                
output_directory = sys.argv[3]               

# Set your NCBI Entrez email  and API key here
Entrez.email = '***REMOVED***'
Entrez.api_key = '92573590fbfb9479ab2167a4f133ee31a408'

def search_genomes(search_term, database='nuccore'):
    """
    Search NCBI GenBank for genome sequences matching the provided search term.

    Args:
    search_term (str): The search term used to query NCBI GenBank.
    database (str): The database to search in (default is 'nuccore').

    Returns:
    list: List of IDs of genome sequences matching the search term.
    """
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
    """
    Fetch GenBank records for the given list of IDs.

    Args:
    id_list (list): List of IDs of genome sequences to fetch.
    database (str): The database to fetch records from (default is 'nuccore').

    Returns:
    list: List of GenBank records fetched for the given IDs.
    """
    handle = Entrez.efetch(db=database, id=",".join(id_list), rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    return list(records)

def save_id_species_mapping(id_species_mapping, mapping_file_path):
    """
    Save the mapping of genome IDs to species names to a file.

    Args:
    id_species_mapping (list): List of tuples containing genome IDs and corresponding species names.
    mapping_file_path (str): Path to the output file where the mapping will be saved.

    Returns:
    None
    """
    with open(mapping_file_path, "w") as f:
        for id, species in id_species_mapping:
            f.write(f"{id}\t{species}\n")

def save_genome_sequences(records, directory, mapping_file_path):
    """
    Save the genome sequences in FASTA format and generate a mapping of genome IDs to species names.

    Args:
    records (list): List of GenBank records containing genome sequences.
    directory (str): Directory where the FASTA files will be saved.
    mapping_file_path (str): Path to the output file where the mapping will be saved.

    Returns:
    None
    """
    id_species_mapping = []
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    for record in records:
        fasta_filepath = os.path.join(directory, f"{record.id}.fasta.gz")
        with gzip.open(fasta_filepath, "wt") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        
        species_name = record.annotations.get('organism', 'Unknown')
        id_species_mapping.append((record.id, species_name))
    
    save_id_species_mapping(id_species_mapping, mapping_file_path)

def main(search_term, output_dir, mapping_file_path):
    """
    Main function to search, fetch, and save genome sequences.

    Args:
    search_term (str): The search term used to query NCBI GenBank.
    output_dir (str): Directory where the genome sequences will be saved.
    mapping_file_path (str): Path to the output file where the mapping will be saved.

    Returns:
    None
    """
    ids = search_genomes(search_term)
    if ids:
        records = fetch_genome_info(ids)
        save_genome_sequences(records, output_dir, mapping_file_path)

if __name__ == "__main__":
    main(query, output_directory, id_species)