from Bio import Entrez

Entrez.email = '***REMOVED***'
Entrez.api_key = '***REMOVED***'

def get_genus_info_from_taxon_id(taxon_id):
    handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
    records = Entrez.read(handle)
    genus_name, genus_taxid = None, None
    # Loop through the LineageEx to find the genus rank and extract name and TaxID
    for lineage_info in records[0]['LineageEx']:
        if lineage_info['Rank'] == 'genus':
            genus_name = lineage_info['ScientificName']
            genus_taxid = lineage_info['TaxId']
            break
    return genus_name, genus_taxid

taxon_id = "81972"  # Example taxon ID for Homo sapiens
genus_name, genus_taxid = get_genus_info_from_taxon_id(taxon_id)
print(f"The genus for taxon ID {taxon_id} is {genus_name} with TaxID {genus_taxid}.")

