import os

def process_files(directory):
    # Path to save the combined file
    output_path = os.path.join(directory, '/home/projects/MAAG/PlantEuka/results/chloroplast/combined_stats_cp.tsv')
    # Start writing to the combined file
    with open(output_path, 'w') as outfile:
        # Write the header once
        outfile.write("Rank\tName\tid\tA\tT\tC\tG\tN\tlength\tA:T %\tC:G %\n")
        
        # Iterate over each file in the directory
        for filename in os.listdir(directory):
            if filename.endswith('.tsv') and 'stats' in filename:
                # Extract Rank and Name from the filename
                parts = filename.split('_')
                rank = parts[0].capitalize()  # Capitalize the first letter
                name = parts[1]
                
                # Build the file path
                file_path = os.path.join(directory, filename)
                
                # Read and process each file
                with open(file_path, 'r') as file:
                    capture = False
                    for line in file:
                        # Start capturing after the id line
                        if line.strip().startswith('id'):
                            capture = True
                            continue
                        # Stop capturing before the Average line
                        if line.strip().startswith('Average'):
                            capture = False
                        # Write the data lines to the output file
                        if capture:
                            outfile.write(f"{rank}\t{name}\t{line}")

# Usage
directory = '/home/projects/MAAG/PlantEuka/results/chloroplast/statistics'
process_files(directory)
