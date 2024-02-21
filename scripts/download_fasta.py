from Bio import Entrez, SeqIO
import os
import gzip

# Set your NCBI Entrez email here
Entrez.email = 'catiacarmobatista@gmail.com'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'

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
    handle = Entrez.efetch(db=database, id=id_list, rettype="fasta", retmode="text")
    records = SeqIO.parse(handle, "fasta")
    return list(records)

def save_genome_sequences(records, directory):
    id_species_mapping = []
    for record in records:
        if not os.path.exists(directory):
            os.makedirs(directory)
        filepath = os.path.join(directory, f"{record.id}.fasta.gz")
        with gzip.open(filepath, "wt") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        # Extract species name from the record description
        species_name = ' '.join(record.description.split(' ')[1:3])  # Adjust as needed
        id_species_mapping.append((record.id, species_name))
    # Call the function to save ID-species mapping
    save_id_species_mapping(id_species_mapping, directory)

def save_id_species_mapping(id_species_mapping, directory):
    filepath = os.path.join(directory, "id_species_mapping.txt")
    with open(filepath, "w") as f:
        for id, species in id_species_mapping:
            f.write(f"{id}\t{species}\n")

def main(search_term, output_dir):
    ids = search_genomes(search_term)
    if ids:
        records = fetch_genome_info(ids)
        save_genome_sequences(records, output_dir)

if __name__ == "__main__":
    search_term = "plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]"
    output_dir = "data/chloroplast"
    main(search_term, output_dir)
