#!/usr/bin/env python3

"""
Description:
This script is designed to download specific gene information from NCBI's GenBank databases.
It fetches genomes based on a user-provided query, extracts the requested gene information, and saves it to a structured directory.
"""

# Libraries 
from Bio import Entrez, SeqIO
import sys
import os

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 2:
    print("Usage: python download_gene.py <organelle>")
    sys.exit(1)

# Extract command-line arguments
organelle = sys.argv[1]
accession_list_file = f"{organelle}/other/accessions.txt"
gene_name = "rbcL"
output_folder = f"genes/{organelle}/{gene_name}"

# Set your NCBI Entrez email and API key here
Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'

def fetch_genome_info(id_list, database='nuccore'):
    """
    Fetch detailed GenBank records for the specified IDs.

    Args:
    id_list (list): A list of IDs for which to retrieve the records.
    database (str): The database from which to fetch records (default is 'nuccore').

    Returns:
    list: List of SeqRecord objects containing genomic information.
    """
    print(f"Fetching genome information for {len(id_list)} IDs...")
    try:
        handle = Entrez.efetch(db=database, id=",".join(id_list), rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "genbank")
        print("Fetch complete.")
        return list(records)
    except Exception as e:
        print(f"Error fetching genome information: {e}")
        return []

def extract_gene_info(records, gene):
    """
    Extract the specified gene information from the provided records.

    Args:
    records (list): List of SeqRecord objects to search within.
    gene (str): Name of the gene to find information about.

    Returns:
    list: A list of tuples containing the ID and location of each found gene.
    """
    print(f"Extracting gene information for gene: {gene}")
    gene_info = []
    for record in records:
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                if gene in feature.qualifiers["gene"]:
                    location = str(feature.location).strip("<>")
                    gene_info.append((record.id, location))
    print(f"Extraction complete. Found {len(gene_info)} gene(s).")
    return gene_info

def save_gene_info(gene_info, output_folder):
    """
    Save the gene information to a structured directory.

    Args:
    gene_info (list): List of tuples containing gene IDs and locations.
    output_folder (str): Path to the directory where the data should be saved.

    """
    print(f"Saving gene information to {output_folder}...")
    try:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        output_file_path = os.path.join(output_folder, "data.txt")
        with open(output_file_path, "w") as f:
            for id, location in gene_info:
                f.write(f"{id} {location}\n")
        print("Save complete.")
    except Exception as e:
        print(f"Error saving gene information: {e}")

def main(accession_list_file, gene_name):
    """
    Main function to manage the workflow of fetching, extracting, and saving gene information.

    Args:
    accession_list_file (str): Path to the file containing accession numbers.
    gene_name (str): Name of the gene to extract information for.
    """
    try:
        print(f"Reading accession list from {accession_list_file}...")
        with open(accession_list_file) as file:
            ids = file.read().splitlines()
        
        if ids:
            print(f"Found {len(ids)} accession ID(s).")
            records = fetch_genome_info(ids)
            gene_info = extract_gene_info(records, gene_name)
            save_gene_info(gene_info, output_folder)
            print(f"{gene_name} gene information saved in {output_folder} directory")
        else:
            print("No accession IDs found.")
    except FileNotFoundError:
        print(f"Error: Accession list file '{accession_list_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main(accession_list_file, gene_name)