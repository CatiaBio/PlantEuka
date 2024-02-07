from Bio import Entrez, SeqIO
import pandas as pd
import os

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'


# Test to get a maximum of 10 sequences 
# def search_genomes(phylum_name, database='nuccore', ret_max=10):
#     search_term = f"{phylum_name}[Orgn] AND chloroplast[Title] AND complete genome[Title]"
#     handle = Entrez.esearch(db=database, term=search_term, retmax=ret_max)
#     record = Entrez.read(handle)
#     handle.close()
#     return record['IdList']

def search_genomes(phylum_name, database='nuccore'):
    search_term = f"{phylum_name}[Orgn] AND chloroplast[Title] AND complete genome[Title]"
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

def extract_taxon_info(genome_records):
    data = []
    for record in genome_records:
        taxon_id = [x for x in record.features if x.type == "source"][0].qualifiers.get("db_xref", ["None"])[0].split(":")[-1]
        data.append({"Accession": record.id, "Organism": record.annotations["organism"], "TaxonID": taxon_id})
    return pd.DataFrame(data)

# To obtain the individual fasta files genomes
# def save_genome_sequences(records, directory="./"):
#     for record in records:
#         SeqIO.write(record, os.path.join(directory, f"{record.id}.fasta"), "fasta")

# To obtain a fasta file with all the genomes 
def save_genome_sequences(records, file_path="all_genomes.fasta"):
    with open(file_path, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

def save_taxon_info(df, filepath="taxon_info.csv"):
    df.to_csv(filepath, index=False)

def main(phylum_name):
    # Create a directory for the output
    output_dir = "sdata/"+f"{phylum_name}_genomes"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Search for genomes
    print(f"Searching for genomes of {phylum_name}...")
    ids = search_genomes(phylum_name)

    # Step 2: Fetch genome information
    print(f"Fetching genome information for {len(ids)} genomes...")
    records = fetch_genome_info(ids)

    # # Step 3: Extract and save taxon information
    # taxon_info_df = extract_taxon_info(records)
    # save_taxon_info(taxon_info_df, os.path.join(output_dir, "taxon_info.csv"))

    # # Step 4: Save genome sequences
    # print(f"Saving genome sequences and taxon information...")
    # save_genome_sequences(records, directory=output_dir)

    # print(f"Done. Data saved to {output_dir}") 

    # Step 3: Extract and save taxon information
    taxon_info_df = extract_taxon_info(records)
    save_taxon_info(taxon_info_df, os.path.join(output_dir, "taxon_info.csv"))

    # Step 4: Save all genome sequences into a single FASTA file
    print("Saving all genome sequences into a single FASTA file...")
    save_genome_sequences(records, file_path=os.path.join(output_dir, "all_genomes.fasta"))

    print(f"Done. Data saved to {output_dir}")
    
if __name__ == "__main__":
    # Example usage: Download chloroplast genomes for the genus "Arabidopsis"
    main("Chlorophyta")
