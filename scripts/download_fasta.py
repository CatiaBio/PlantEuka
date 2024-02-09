from Bio import Entrez, SeqIO
import os

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

def search_genomes(genus_name, database='nuccore'):
    search_term = f"{genus_name}[Orgn] AND chloroplast[Title] AND complete genome[Title]"
    id_list = []
    
    # Initial search to get count
    handle = Entrez.esearch(db=database, term=search_term, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    count = int(record['Count'])
    
    batch_size = 500  # NCBI recommends not to exceed 500 for large queries
    for start in range(0, count, batch_size):
        handle = Entrez.esearch(db=database, term=search_term, retstart=start, retmax=batch_size)
        record = Entrez.read(handle)
        handle.close()
        id_list.extend(record['IdList'])
    
    return id_list

def fetch_genome_info(id_list, database='nuccore'):
    handle = Entrez.efetch(db=database, id=id_list, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "gb")
    return list(records)

def save_genome_sequences(records, directory="./"):
    for record in records:
        SeqIO.write(record, os.path.join(directory, f"{record.id}.fasta"), "fasta")

def read_genus_names_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def main(file_path):
    genus_names = read_genus_names_from_file(file_path)
    for genus_name in genus_names:
        output_dir = "data/"+f"{genus_name}_genomes"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        print(f"Searching for genomes of {genus_name}...")
        ids = search_genomes(genus_name)
        print(f"Fetching genome information for {len(ids)} genomes...")
        records = fetch_genome_info(ids)
        print(f"Saving genome sequences and taxon information...")
        save_genome_sequences(records, directory=output_dir)
        print(f"Done. Data saved to {output_dir}")

if __name__ == "__main__":
    file_path = 'data/plant_genera_taxids_names.tsv'  # Replace with your actual text file path
    main(file_path)