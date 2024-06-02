#!/usr/bin/env python3

"""
Description:
This script is designed to download specific gene information from NCBI's GenBank databases.
It fetches genomes based on a user-provided query, extracts the requested gene information, and saves it to a structured directory.

Usage:
./download_gene.py <accession list file> <gene name> <gene location>

Arguments:
<accession list file>: Path to the file containing accession numbers.
<gene name>: Name of the gene to extract information for (e.g., rbcL, psaA).
<organelle type>: Location to categorize and save the gene data (e.g., chloroplast).

Example usage following PlantEuka folder organization:
./scripts/download_gene.py missing_accessions.txt rbcL chloroplast
"""

# Libraries 
from Bio import Entrez, SeqIO
import sys
import os

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 4:
    print("Usage: ./download_gene.py <accession list file> <gene name> <organelle type>")
    sys.exit(1)

# Extract command-line arguments
accession_list_file = sys.argv[1]
gene_name = sys.argv[2]
organelle_type = sys.argv[3]
output_folder = f"genes/{organelle_type}/{gene_name}"

# Set your NCBI Entrez email and API key here
Entrez.email = 'catiacarmobatista@gmail.com'
Entrez.api_key = '92573590fbfb9479ab2167a4f133ee31a408'

def fetch_genome_info(id_list, database='nuccore'):
    """
    Fetch detailed GenBank records for the specified IDs.

    Args:
    id_list (list): A list of IDs for which to retrieve the records.
    database (str): The database from which to fetch records (default is 'nuccore').

    Returns:
    list: List of SeqRecord objects containing genomic information.
    """
    handle = Entrez.efetch(db=database, id=",".join(id_list), rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    return list(records)

def extract_gene_info(records, gene):
    """
    Extract the specified gene information from the provided records.

    Args:
    records (list): List of SeqRecord objects to search within.
    gene (str): Name of the gene to find information about.

    Returns:
    list: A list of tuples containing the ID and location of each found gene.
    """
    gene_info = []
    for record in records:
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                if gene in feature.qualifiers["gene"]:
                    location = str(feature.location).strip("<>")
                    gene_info.append((record.id, location))
    return gene_info

def save_gene_info(gene_info, output_folder):
    """
    Save the gene information to a structured directory.

    Args:
    gene_info (list): List of tuples containing gene IDs and locations.
    output_folder (str): Path to the directory where the data should be saved.

    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_file_path = os.path.join(output_folder, "data.txt")
    with open(output_file_path, "w") as f:
        for id, location in gene_info:
            f.write(f"{id} {location}\n")

def main(accession_list_file, gene_name, gene_location):
    """
    Main function to manage the workflow of fetching, extracting, and saving gene information.

    Args:
    accession_list_file (str): Path to the file containing accession numbers.
    gene_name (str): Name of the gene to extract information for.
    gene_location (str): Location to categorize and save the gene data.
    """
    with open(accession_list_file) as file:
        ids = file.read().splitlines()
    
    if ids:
        records = fetch_genome_info(ids)
        gene_info = extract_gene_info(records, gene_name)
        save_gene_info(gene_info, f"genes/{gene_location}/{gene_name}")
        print(f"{gene_name} gene information saved in {gene_location} directory")

if __name__ == "__main__":
    main(accession_list_file, gene_name, organelle_type)
