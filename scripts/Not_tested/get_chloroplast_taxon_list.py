from Bio import Entrez, SeqIO
import pandas as pd
import os

# Set your NCBI Entrez email here
Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'

def search_genomes(database='nuccore', ret_max=500):
    search_term = "chloroplast[Title] AND complete genome[Title] AND (\"plants\"[Filter]) AND refseq[Filter]"
    handle = Entrez.esearch(db=database, term=search_term, retmax=ret_max)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def fetch_taxon_info(id_list, database='nuccore'):
    data = []
    for batch in chunks(id_list, 500):  # Process in batches
        handle = Entrez.efetch(db=database, id=','.join(batch), rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        for record in records:
            # Safely extract taxonomy ID, organism, and accession
            try:
                # For XML format, adjust paths according to actual structure
                accession = record['GBSeq_locus']  # Example, adjust as needed
                organism = record['GBSeq_organism']
                # Taxonomy ID might need fetching from another field or additional parsing
                taxon_id = 'N/A'  # Placeholder, adjust your approach to extract it
                for feature in record['GBSeq_feature-table']:
                    if feature['GBFeature_key'] == "source":
                        for qualifier in feature['GBFeature_quals']:
                            if qualifier['GBQualifier_name'] == "db_xref":
                                taxon_value = qualifier['GBQualifier_value']
                                if taxon_value.startswith("taxon:"):
                                    taxon_id = taxon_value.split(":")[1]
                data.append({"Accession": accession, "Organism": organism, "TaxonID": taxon_id})
            except KeyError as e:
                print(f"Key error: {e} in record {record['GBSeq_locus']}")
        handle.close()
    return pd.DataFrame(data)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def save_taxon_info(df, filepath="taxon_info.tsv"):
    df.to_csv(filepath, sep='\t', index=False)

def main():
    # Step 1: Search for chloroplast complete genomes
    print("Searching for chloroplast complete genomes...")
    ids = search_genomes(ret_max=20)  # Request up to 500 ids
    # ids = search_genomes()

    # Step 2: Fetch and save taxon information for the first 500 records
    print(f"Fetching taxon information for {len(ids)} records...")
    taxon_info_df = fetch_taxon_info(ids[:500])  # Ensure only the first 500 are processed

    # Save to TSV
    output_file = "data/chloroplast_taxon_info.tsv"
    print(f"Saving taxon information to {output_file}...")
    save_taxon_info(taxon_info_df, filepath=output_file)

    print("Done.")

if __name__ == "__main__":
    main()
