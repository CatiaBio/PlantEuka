from Bio import Entrez, SeqIO
import csv

# Set your NCBI Entrez email and API key
Entrez.email = 'catiacarmobatista@gmail.com'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def get_id_list(file_path):
    id_list = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            id_list.append(row[0])  # Fixed to use append for list
    return id_list

def search_and_download_genomes(id_list, database='nucleotide', output_file='gene_sequences.fasta'):
    sequences = []
    for id in id_list:
        search_term = f"{id}[Accession] AND rbcL[Gene]"  
        handle = Entrez.esearch(db=database, term=search_term, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()
        
        # Proceed if there are search results
        if search_results['IdList']:
            seq_id = search_results['IdList'][0]  # Assuming we take the first match
            fetch_handle = Entrez.efetch(db=database, id=seq_id, rettype="fasta", retmode="text")
            record = SeqIO.read(fetch_handle, "fasta")
            fetch_handle.close()
            sequences.append(record)
    
    # Save all sequences into a single file
    with open(output_file, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    print(f"All sequences have been saved to {output_file}")

if __name__ == "__main__":
    file_path = '../data/mapping_id_taxID_cp.txt'
    id_list = get_id_list(file_path)
    search_and_download_genomes(id_list)
