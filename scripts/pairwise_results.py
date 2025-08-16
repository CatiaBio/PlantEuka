
#!/usr/bin/env python3

import os
import gzip
import re

# For Snakemake compatibility
directory = snakemake.input.results_dir
output_file = snakemake.output.results_file

# Header for the TSV file
header = 'rank\tname\tpair1\tpair2\tlength\tidentity\tsimilarity\tgaps\tscore\n'

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Open the output file and write the header
with open(output_file, 'w') as outfile:
    outfile.write(header)

    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.stretcher.gz'):  # Check for compressed needle files
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
                outfile.write(f"{rank}\t{name}\t{pair1}\t{pair2}\t{data.get('Length', 'N/A')}\t{data.get('Identity', 'N/A')}\t{data.get('Similarity', 'N/A')}\t{data.get('Gaps', 'N/A')}\t{data.get('Score', 'N/A')}\n")

print("Data extraction complete.")