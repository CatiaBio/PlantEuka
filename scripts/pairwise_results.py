#!/usr/bin/env python3

import os
import gzip
import re
import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 2:
    print("Usage: python3 download_genomes.py <organelle>")
    sys.exit(1)

# Extract command-line arguments
organelle = sys.argv[1]

# Directory containing the files
directory = f'{organelle}/results/pairwise/stretcher'

# Output TSV file
output_file = '{organelle}/results/pairwise/pairwise_results.tsv'

# Header for the TSV file
header = 'Rank\tName\tpair1\tpair2\tLength\tIdentity\tSimilarity\tGaps\tScore\n'

# Open the output file and write the header
with open(output_file, 'w') as outfile:
    outfile.write(header)

    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.needle.gz'):  # Check for compressed needle files
            filepath = os.path.join(directory, filename)
            # Extract the rank and name from the filename
            parts = filename.replace('.needle.gz', '').split('_')
            rank = parts[0]  # Assuming the rank is the first part
            name = parts[1]  # Assuming the name is the second part

            with gzip.open(filepath, 'rt') as file:  # Open the gzip file in text mode
                data = {}  # Dictionary to store the data
                capture_pairs = False  # Flag to start capturing pair names
                
                # Read through each line in the file
                for line in file:
                    if 'Aligned_sequences: 2' in line:
                        capture_pairs = True
                    elif capture_pairs and line.startswith('# 1:'):
                        pair1 = line.split(':')[1].strip()
                    elif capture_pairs and line.startswith('# 2:'):
                        pair2 = line.split(':')[1].strip()
                        capture_pairs = False  # Stop capturing after finding both pairs
                    elif 'Length:' in line:
                        data['Length'] = line.split(':')[-1].strip()
                    elif 'Identity:' in line:
                        identity_results = re.findall(r"\((\d+\.\d+)%\)", line)
                        data['Identity'] = identity_results[0] if identity_results else 'N/A'
                    elif 'Similarity:' in line:
                        similarity_results = re.findall(r"\((\d+\.\d+)%\)", line)
                        data['Similarity'] = similarity_results[0] if similarity_results else 'N/A'
                    elif 'Gaps:' in line:
                        gaps_results = re.findall(r"\(\s*(\d+\.\d+)%\)", line)
                        data['Gaps'] = gaps_results[0] if gaps_results else 'N/A'
                    elif 'Score:' in line:
                        data['Score'] = line.split(':')[-1].strip()

                # Write the extracted data to the TSV file
                outfile.write(f"{rank}\t{name}\t{pair1}\t{pair2}\t{data.get('length', 'N/A')}\t{data.get('identity', 'N/A')}\t{data.get('similarity', 'N/A')}\t{data.get('gaps', 'N/A')}\t{data.get('score', 'N/A')}\n")
                
print("Data extraction complete.")