#!/bin/bash


   # Example of how you might pass the pairs to another script or command
        # ./process_pair.sh "$pair1" "$pair2"
        # or simply:
        # echo "$pair1 $pair2" | parallel --some-options "command {1} {2}"
        # ./process_pair.sh "$pair1" "$pair2" | parallel | parallel --slf list_serv_2_14 "nice -19 /home/ctools/EMBOSS-6.6.0/emboss/stretcher -asequence {pair1} -bsequence {pair2} -gapopen 16 -gapextend 4 -outfile {pair1_pair2}.needle" &


# Define the directory where your list files are stored
list_dir="$(pwd)/other/lists"

# Ensure the list directory exists
if [ ! -d "$list_dir" ]; then
    echo "List directory does not exist: $list_dir"
    exit 1
fi

# Loop through all list files in the specified directory
for list_file in "${list_dir}"/*_pairs.txt; do
    echo "Processing list file: $list_file"

    # Read each line in the list file
    while IFS= read -r line; do
        # Split the line into pair1 and pair2
        pair1=$(echo "$line" | cut -d ' ' -f1)
        pair2=$(echo "$line" | cut -d ' ' -f2)

        # Output the pairs or pass them to another command
        echo "Pair1: $pair1"
        echo "Pair2: $pair2"

    done < "$list_file"
done

echo "All list files have been processed."
