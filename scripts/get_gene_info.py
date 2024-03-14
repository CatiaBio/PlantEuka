from Bio import Entrez
import csv

# Set your NCBI Entrez email and API key
Entrez.email = '***REMOVED***'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def download_gene_data(accession_number, gene_name, output_file):
    """
    Download tabular data for a specific gene from NCBI for a given accession number.

    Parameters:
    accession_number (str): The accession number of interest.
    gene_name (str): The name of the gene of interest (e.g., "psaA").
    output_file (str): Path to the output file where the data will be saved.
    """
    # Formulate the query
    query = f"{gene_name}[Gene Name] AND {accession_number}[Accession]"
    
    # Search the NCBI gene database
    search_results = Entrez.read(Entrez.esearch(db="gene", term=query, retmax=10))
    id_list = search_results['IdList']
    
    if id_list:
        # Fetch the detailed gene information based on the search results
        fetch_results = Entrez.efetch(db="gene", id=",".join(id_list), retmode="text", rettype="tabular")
        data = fetch_results.read()
        
        # For simplicity, this example just writes the raw XML data to the output file.
        # You might want to parse the XML and extract specific information instead.
        # Inside your download_gene_data function, before writing data to the file:
        if isinstance(data, bytes):
            data = data.decode('utf-8')  # Decode bytes to str
        with open(output_file, 'w') as file:
            file.write(data)

        print(f"Data for {gene_name} gene and accession {accession_number} saved to {output_file}")
    else:
        print(f"No data found for {gene_name} gene with accession {accession_number}")

if __name__ == "__main__":
    # Path to your TSV file containing accession numbers and species names
    tsv_file_path = '../data/mapping_id_species_cp.txt'  
    
    # Specify the gene of interest
    gene_name = "psaA"
    
    with open(tsv_file_path) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            accession_number, species_name = row
            # Define an output file for each accession number
            output_file = f"{species_name.replace(' ', '_')}_{accession_number}_{gene_name}_data.txt"
            download_gene_data(accession_number, gene_name, output_file)
