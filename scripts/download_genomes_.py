import os
import time
from Bio import Entrez
import gzip
import requests

# Access paths from Snakemake
accession_list_file = snakemake.input.accession_list  # List of all accessions
downloaded_list_file = snakemake.output.updated_downloaded_list  # Track already downloaded accessions
ncbi_info_file = snakemake.input.ncbi_info  # NCBI info file for email and API key
genomes_dir = snakemake.params.genomes_dir  # Directory to store downloaded genome files

# Read email and API key from ncbi_info.txt
with open(ncbi_info_file, "r") as info:
    get_info = info.readlines()
    email = get_info[0].strip()
    api_key = get_info[1].strip()

# Set Entrez email and API key
Entrez.email = email
Entrez.api_key = api_key

# Ensure the genomes directory exists
if not os.path.exists(genomes_dir):
    os.makedirs(genomes_dir)

# Read accession numbers from accession_list.txt
with open(accession_list_file, "r") as file:
    all_accessions = set(line.strip() for line in file)

# Check for downloaded accessions if the file exists
if os.path.exists(downloaded_list_file):
    with open(downloaded_list_file, "r") as file:
        downloaded_accessions = set(line.strip() for line in file)
else:
    downloaded_accessions = set()

# Identify accessions that need to be downloaded
missing_accessions = all_accessions - downloaded_accessions

# Function to download and save a genome as a .fasta.gz file with API key support
def download_genome(accession):
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={accession}&db=nuccore&report=fasta&extrafeat=null&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=10000000&api_key={api_key}"
    
    # Delay to respect NCBI's rate limiting if required
    time.sleep(0.3)  # 300ms between requests

    response = requests.get(url)
    if response.status_code == 200:
        # Save the response as a .fasta.gz file
        fasta_path = os.path.join(genomes_dir, f"{accession}.fasta.gz")
        with gzip.open(fasta_path, "wt") as f:
            f.write(response.text)
        return True
    else:
        print(f"Failed to download {accession}: HTTP {response.status_code}")
        return False

# Download missing genomes and update the downloaded_accessions_list.txt file
newly_downloaded = []
for accession in missing_accessions:
    if download_genome(accession):
        newly_downloaded.append(accession)

# Update the downloaded_accessions_list.txt with newly downloaded accessions
if newly_downloaded:
    with open(downloaded_list_file, "a") as file:
        for accession in newly_downloaded:
            file.write(accession + "\n")
    print(f"Downloaded {len(newly_downloaded)} new genome(s) and updated {downloaded_list_file}.")
else:
    print("No new genomes to download.")