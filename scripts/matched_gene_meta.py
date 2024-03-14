import csv

# Paths to your TSV files
file1_path = '../data/rubisco_gene_metadata.txt'  # Replace with the path to your first TSV file
file2_path = '../data/mapping_id_species_cp.txt'  # Replace with the path to your second TSV file
output_file_path = '../data/matched_records.tsv'  # Output file to save the matching records

# Read the accessions from file 2 into a set for fast lookup
accessions_file2 = set()
with open(file2_path, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        accessions_file2.add(row[0])  # Assuming the accession is the first element

# Now read file 1 and write out the matching records to the output file
with open(file1_path, 'r') as file1, open(output_file_path, 'w', newline='') as output_file:
    reader = csv.reader(file1, delimiter='\t')
    writer = csv.writer(output_file, delimiter='\t')
    
    #Optional: Write header to the output file if needed
    header = next(reader)  # Uncomment if your file has a header
    writer.writerow(header)  # Uncomment to write the header to the output file
    
    for row in reader:
        # Check if the 'genomic_nucleotide_accession.version' column matches any accession from file 2
        if row[11] in accessions_file2:  # Adjust the index 8 if necessary to match 'genomic_nucleotide_accession.version' column
            writer.writerow(row)

print("Matching records have been saved to:", output_file_path)
