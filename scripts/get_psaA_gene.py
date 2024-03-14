from Bio import Entrez, SeqIO
import time

# Set your NCBI Entrez email and API key
Entrez.email = 'catiacarmobatista@gmail.com'
Entrez.api_key = '03d4de8e8f44ed94183a4ba257fab9752709'


def download_psaA_sequences(species_name, output_filename_prefix):
    """
    Download psaA gene sequences for a given species from NCBI and save as FASTA.

    :param species_name: The common or scientific name of the species.
    :param output_filename_prefix: Prefix for the output file names.
    """
    query = f'psaA[Gene Name] AND "{species_name}"[Organism] AND chloroplast[Filter]'
    search_handle = Entrez.esearch(db='nucleotide', term=query, retmax=10)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    id_list = search_results['IdList']
    if id_list:
        fetch_handle = Entrez.efetch(db='nucleotide', id=id_list, rettype='fasta', retmode='text')
        records = list(SeqIO.parse(fetch_handle, 'fasta'))
        fetch_handle.close()

        if records:
            output_filename = f'{output_filename_prefix}_{species_name.replace(" ", "_")}_psaA.fasta'
            SeqIO.write(records, output_filename, 'fasta')
            print(f'Saved: {output_filename}')
        else:
            print(f'No records found for {species_name}')
    else:
        print(f'No results found for {species_name}')

    time.sleep(1)  # to avoid overwhelming NCBI servers

# Read species list from TSV file
tsv_file_path = '../data/mapping_id_species_cp.txt'  # Update this to your TSV file path
with open(tsv_file_path, 'r') as file:
    for line in file:
        accession, species = line.strip().split('\t')
        download_psaA_sequences(species, accession)

print("Download complete.")
