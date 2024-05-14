import subprocess

def extract_taxid(accession_file, db_file, output_file):
    """
    Extracts taxids for accession versions from a gzipped NCBI accession2taxid file.
    
    Args:
    accession_file (str): Path to the file containing accession.version numbers.
    db_file (str): Path to the gzipped nucl_gb.accession2taxid.gz file.
    output_file (str): Path to save the output list of accession.version and taxid.
    """
    # Open the output file for writing
    with open(output_file, 'w') as outfile:
        outfile.write("accession.version\ttaxid\n")  # Writing header to the output file
        
        # Read accession versions from the input file
        with open(accession_file, 'r') as infile:
            for accession in infile:
                accession = accession.strip()
                if accession:
                    # Use zgrep to search the gzipped database file
                    result = subprocess.run(['zgrep', f"^{accession}\t", db_file],
                                            capture_output=True, text=True)
                    
                    # Process the output from zgrep
                    if result.stdout:
                        # Extract the taxid from the output line (expected format: accession\taccession.version\ttaxid\tgi)
                        fields = result.stdout.split()
                        taxid = fields[2] if len(fields) > 2 else 'Not_Found'
                        outfile.write(f"{accession}\t{taxid}\n")
                    else:
                        # No result found, write accession with 'Not_Found' taxid
                        outfile.write(f"{accession}\tNot_Found\n")

# Usage example
if __name__ == "__main__":
    extract_taxid('other/2405_all_cp_genomes.txt', 'nucl_gb.accession2taxid.gz', '2405_id_taxid_cp.txt')