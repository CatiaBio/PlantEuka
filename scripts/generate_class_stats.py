#!/usr/bin/env python3

import csv

def load_lineage_data(lineage_file):
    """
    Loads lineage data to create a detailed map of ID to all lineage information.
    """
    lineage_dict = {}
    with open(lineage_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            lineage_dict[row['ID']] = row
    return lineage_dict

def process_combined_stats(combined_stats_file, lineage_dict, output_file):
    """
    Processes the combined stats file and writes the output with full lineage information plus selected genomic stats.
    """
    with open(combined_stats_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        csv_reader = csv.DictReader(infile, delimiter='\t')
        first_row = next(csv_reader)  # Read the first row to get additional headers
        stats_headers = ['A', 'T', 'C', 'G', 'N', 'length', 'A:T %', 'C:G %']

        # Define new headers combining lineage and stats
        lineage_headers = ['ID', 'TaxID', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum']
        all_headers = lineage_headers + stats_headers
        csv_writer = csv.DictWriter(outfile, fieldnames=all_headers, delimiter='\t')
        csv_writer.writeheader()

        # Reset the reader to the start of the file and process all rows
        infile.seek(0)
        next(csv_reader)  # Skip header row

        for row in csv_reader:
            id = row['id']
            if id in lineage_dict:
                # Start with copying all lineage details
                output_row = {**lineage_dict[id]}
                # Add selected genomic stats
                for header in stats_headers:
                    output_row[header] = row[header]
                csv_writer.writerow(output_row)

if __name__ == "__main__":
    # Define file paths
    lineage_file = 'other/lineage_with_id.tsv'
    combined_stats_file = 'results/combined_stats_cp.tsv'
    output_file = 'results/combined_class_stats_cp.tsv'

    # Load lineage data
    lineage_dict = load_lineage_data(lineage_file)

    # Process combined stats and output results
    process_combined_stats(combined_stats_file, lineage_dict, output_file)
    print("Output has been written to", output_file)
