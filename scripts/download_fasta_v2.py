from Bio import Entrez, SeqIO
import os
import gzip
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Download fasta files.')
parser.add_argument('--search_str', required=True, help='Query string')
parser.add_argument('--mapping_file', required=True, help='Path to save the mapping file')
parser.add_argument('--output_dir', required=True, help='Path to the output directory')
args = parser.parse_args()

# Use argparse to parse command line arguments
search_term = args.search_str
mapping_file_path = args.mapping_file
output_dir = args.output_dir

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'

def search_genomes(search_term, database='nuccore'):
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
    # Fetch GenBank records for the given list of IDs
    handle = Entrez.efetch(db=database, id=",".join(id_list), rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    return list(records)

def save_id_species_mapping(id_species_mapping, mapping_file_path):
    with open(mapping_file_path, "w") as f:
        for id, species in id_species_mapping:
            f.write(f"{id}\t{species}\n")

def save_genome_sequences(records, directory, mapping_file_path):
    id_species_mapping = []
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    for record in records:
        # fasta_filepath = os.path.join(directory, f"{record.id}.fasta.gz")
        # with gzip.open(fasta_filepath, "wt") as output_handle:
        #     SeqIO.write(record, output_handle, "fasta")
        
        species_name = record.annotations.get('organism', 'Unknown')
        id_species_mapping.append((record.id, species_name))
    
    save_id_species_mapping(id_species_mapping, mapping_file_path)



def main(search_term, output_dir, mapping_file_path):
    ids = search_genomes(search_term)
    if ids:
        records = fetch_genome_info(ids)
        save_genome_sequences(records, output_dir, mapping_file_path)

if __name__ == "__main__":
    main(search_term, output_dir, mapping_file_path)
