from Bio import Entrez
import csv
import os

# Set your NCBI Entrez email and API key
Entrez.email = 'catiacarmobatista@gmail.com'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def initialize_output_file_with_header(output_file):
    """
    Initializes the output file with an appropriate header by performing a dummy fetch operation.
    Assumes all fetched entries have a consistent header format.
    """
    # Perform a dummy fetch to get the header; you may replace 'NC_000852' with a known, consistent ID
    dummy_fetch_result = Entrez.efetch(db="gene", id="NC_000852", retmode="text", rettype="tabular")
    dummy_data = dummy_fetch_result.read()
    if isinstance(dummy_data, bytes):
        dummy_data = dummy_data.decode('utf-8')
    
    # Write the header line to the output file
    with open(output_file, 'w') as file:
        header_line = dummy_data.split('\n')[0] + '\n'  # Extract and write only the header
        file.write(header_line)

def download_gene_data(accession_number, gene_name, output_file, has_gene_file, no_gene_file):
    """
    Download tabular data for a specific gene from NCBI for a given accession number
    and append it to a single output file. Log accession numbers to 'has_gene.tsv'
    or 'no_gene.tsv' based on whether data was found.
    """
    query = f"{gene_name}[Gene Name] AND {accession_number}[Accession]"
    search_results = Entrez.read(Entrez.esearch(db="gene", term=query, retmax=10))
    id_list = search_results['IdList']
    
    if id_list:
        fetch_results = Entrez.efetch(db="gene", id=",".join(id_list), retmode="text", rettype="tabular")
        data = fetch_results.read()
        
        if isinstance(data, bytes):
            data = data.decode('utf-8')
        
        # Append the data without the header
        data_lines = data.split('\n')[1:]  # Skip the header line
        with open(output_file, 'a') as file:
            file.write('\n'.join(data_lines))
            file.write('\n')  # Ensure separation between entries
        
        with open(has_gene_file, 'a') as file:
            file.write(f"{accession_number}\n")
    else:
        with open(no_gene_file, 'a') as file:
            file.write(f"{accession_number}\n")

if __name__ == "__main__":
    output_dir = 'genes/chloroplast/rbcL'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file = os.path.join(output_dir, "all_genes_data.tsv")
    has_gene_file = os.path.join(output_dir, "has_gene.tsv")
    no_gene_file = os.path.join(output_dir, "no_gene.tsv")

    # Initialize the output file with the header before appending data
    initialize_output_file_with_header(output_file)

if __name__ == "__main__":
    # Path to your TSV file containing accession numbers and species names
    tsv_file_path = 'genomes/chloroplast/id_taxid_cp.tsv'
    
    # Output directory
    output_dir = 'genes/chloroplast/rbcL'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "all_genes_data.tsv")
    has_gene_file = os.path.join(output_dir, "has_gene.tsv")
    no_gene_file = os.path.join(output_dir, "no_gene.tsv")

    # Initialize the output file with the header before appending data
    initialize_output_file_with_header(output_file)
    
    # Specify the gene of interest
    gene_name = "rbcL"
    
    # Open the log files once to ensure they're empty before appending data
    open(has_gene_file, 'w').close()
    open(no_gene_file, 'w').close()
    
    with open(tsv_file_path) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            accession_number, species_name = row
            # Append data for each accession number to the single output file and log IDs
            download_gene_data(accession_number, gene_name, output_file, has_gene_file, no_gene_file)