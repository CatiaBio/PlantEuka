#!/bin/bash

# Define the directory where your list files are stored
list_dir="$(pwd)/other/lists/test"

# Define the directories for the outputs and logs
results_dir="$(pwd)/results/pairwise"
logs_dir="$(pwd)/logs"

# Ensure the list, results, and logs directories exist
mkdir -p "$list_dir" "$results_dir" "$logs_dir" 

if [ ! -d "$list_dir" ]; then
    echo "List directory does not exist: $list_dir"
    exit 1
fi

# Loop through all list files in the specified directory
for list_file in "${list_dir}"/*_pairs.txt; do
    echo "Processing list file: $list_file"
    
    # Define output file for the parallel input
    parallel_file="${list_file%_pairs.txt}_parallel.txt"
    
    # Open file descriptor to write parallel commands
    exec 3>"$parallel_file"

    # Read each line from the list file and modify it to generate the necessary parallel input
    while IFS= read -r line; do
        # Replace '.fasta.gz' with '.fasta' in each path
        pair1=$(echo "$line" | cut -d ' ' -f1 | sed 's/\.gz$//')
        pair2=$(echo "$line" | cut -d ' ' -f2 | sed 's/\.gz$//')
        
        # Construct the taxonomy and genus prefix from the directory structure
        path_entry=(${line//\// }) 
        # For absolute paths in the server
        taxonomy_genus=$(echo ${path_entry[7]})
        taxonomy_species=$(echo ${path_entry[8]})

        # For path my computer 
        #taxonomy_genus=$(echo ${path_entry[3]})
        #taxonomy_species=$(echo ${path_entry[4]})

        pair1_name=$(basename "${line%% *}" .fasta.gz)
        pair2_name=$(basename "${line##* }" .fasta.gz)

        # Define output and log file names
        output_name="${taxonomy_genus}_${taxonomy_species}_${pair1_name}_${pair2_name}.needle"
        log_name="${taxonomy_genus}_${taxonomy_species}_${pair1_name}_${pair2_name}.log"
        
        # Write the modified line with full paths to the results and log directories
        echo "$pair1 $pair2 $results_dir/$output_name $logs_dir/$log_name" >&3
    
    done < "$list_file"

    # Close the file descriptor
    exec 3>&-

    echo "Parallel commands written to $parallel_file"
done