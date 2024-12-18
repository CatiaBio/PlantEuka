#!/usr/bin/env python3

"""
Description:
This script is designed to download genome sequences from the NCBI database based on a provided search query. 
It retrieves genome records matching the specified search query, saves the sequences in FASTA format, 
and generates a file listing the accession numbers of the downloaded genomes.

Usage:
python3 download_genomes.py 
"""

# Libraries
from Bio import Entrez, SeqIO
from snakemake import snakemake
import os
import sys
import gzip

# Access the file paths from Snakemake's input and output
input_file = snakemake.input.config  
output_file = snakemake.output[0]    

with open(ncbi_info, "r") as info:
        get_info = info.readlines()
        email = get_info[0].strip()
        api_key = get_info[1].strip()

# Extract command-line arguments
accession_version_file = output_file
organelle = "chloroplast"
query = f"plants[filter] AND refseq[filter] AND {organelle}[filter] AND complete genome[Title]"                 
accession_version_file = f"{organelle}/{output_file}"                
output_directory = f"{organelle}/genomes/original"              

# Set your NCBI Entrez email and API key here
Entrez.email = email
Entrez.api_key = api_key

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
    handle = Entrez.esearch(db=database, term=search_term, retmax=1)
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

def save_accession_numbers(records, mapping_file_path):
    """
    Save the accession numbers to a file.

    Args:
    records (list): List of GenBank records containing genome sequences.
    mapping_file_path (str): Path to the output file where the accession numbers will be saved.

    Returns:
    None
    """
    directory = os.path.dirname(mapping_file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    with open(mapping_file_path, "w") as f:
        for record in records:
            f.write(f"{record.id}\n")

def save_genome_sequences(records, directory, mapping_file_path):
    """
    Save the genome sequences in FASTA format and generate a file with the accession numbers.

    Args:
    records (list): List of GenBank records containing genome sequences.
    directory (str): Directory where the FASTA files will be saved.
    mapping_file_path (str): Path to the output file where the accession numbers will be saved.

    Returns:
    None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    for record in records:
        fasta_filepath = os.path.join(directory, f"{record.id}.fasta.gz")
        with gzip.open(fasta_filepath, "wt") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
    
    save_accession_numbers(records, mapping_file_path)

def main(search_term, output_dir, mapping_file_path):
    """
    Main function to search, fetch, and save genome sequences.

    Args:
    search_term (str): The search term used to query NCBI GenBank.
    output_dir (str): Directory where the genome sequences will be saved.
    mapping_file_path (str): Path to the output file where the accession numbers will be saved.

    Returns:
    None
    """
    ids = search_genomes(search_term)
    if ids:
        records = fetch_genome_info(ids)
        save_genome_sequences(records, output_dir, mapping_file_path)

if __name__ == "__main__":
    main(query, output_directory, accession_version_file)
