#!/bin/bash

# Example of how you might pass the pairs to another script or command
# nohup ./scripts/run_pairs_parallel.sh | parallel --slf list_serv_2_14 --colsep ' ' "nice -19 /home/ctools/EMBOSS-6.6.0/emboss/stretcher -asequence {1} -bsequence {2} -gapopen 16 -gapextend 4 -outfile {3}" > {4} 2>&1 &


# /genomes/chloroplast/sorted/genus/Abies/NC_045884.1.fasta
# /home/projects/MAAG/PlantEuka/genomes/chloroplast/sorted/genus/Abies/NC_026892.1.fasta.gz

# Define the directory where your list files are stored
list_dir="$(pwd)/other/lists/combined"

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
    
    # Read each line from the list file and modify it to generate the necessary parallel input
    while IFS= read -r line; do
        
        # Extract the base names without the extension and path to generate output and log filenames
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
        
        # Output the modified line with the paths without .gz and the new file names
        echo "$pair1 $pair2 $results_dir/$output_name $logs_dir/$log_name"
    
    done < "$list_file"
done