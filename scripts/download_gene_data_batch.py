from Bio import Entrez
import csv
import os
import time

# Set your NCBI Entrez email and API key
Entrez.email = '***REMOVED***'  # Use your actual email
Entrez.api_key = '***REMOVED***'  # Use your actual API key

def download_gene_data_batch(taxids, gene_name, output_file):
    """
    Download gene data for a list of taxonomic identifiers (taxids) in batches
    and append it to a single output file, with a pause between each batch to respect rate limits.

    Parameters:
    taxids (list): The list of taxonomic identifiers of interest.
    gene_name (str): The name of the gene of interest (e.g., "rbcL").
    output_file (str): Path to the output file where the data will be appended.
    """
    batch_size = 500
    first_entry = True  # Flag to check if it's the first entry across all batches

    for start in range(0, len(taxids), batch_size):
        batch = taxids[start:start+batch_size]
        
        for taxid in batch:
            query = f"{gene_name}[Gene Name] AND {taxid}[Organism:exp]"
            search_results = Entrez.read(Entrez.esearch(db="gene", term=query, retmax=1))
            id_list = search_results['IdList']

            if id_list:
                fetch_results = Entrez.efetch(db="gene", id=",".join(id_list), retmode="text", rettype="tabular")
                data = fetch_results.read()

                if isinstance(data, bytes):
                    data = data.decode('utf-8')

                data_lines = data.split('\n')
                with open(output_file, 'a') as file:
                    if first_entry:  # Include header for the first entry
                        file.write(data + '\n')
                        first_entry = False
                    else:  # Exclude the header for subsequent entries
                        file.write('\n'.join(data_lines[1:]))
                        file.write('\n')

        print(f"Completed a batch ending with taxid {taxid}")
        if start + batch_size < len(taxids):  # Check if there are more batches
            print("Pausing for 10 seconds before the next batch...")
            time.sleep(10)  # Pause for 10 seconds

if __name__ == "__main__":
    output_dir = 'genes/chloroplast/rbcL'
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "all_genes_data.tsv")
    gene_name = "rbcL"

    # Ensure the output file is empty before starting
    open(output_file, 'w').close()

    # Read taxids and species names from the TSV file
    taxids = []
    with open('genomes/chloroplast/id_taxid_cp.tsv') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            _, taxid = row  # Assuming each row is "accession_number\ttaxid"
            taxids.append(taxid)

    # Download data for taxids in batches with a pause between each
    download_gene_data_batch(taxids, gene_name, output_file)
