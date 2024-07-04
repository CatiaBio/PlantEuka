#!/bin/bash

# Get the current working directory
base_dir=$(pwd)

# Path to the list of files
file_list="${base_dir}/chloroplast/other/vg_paths.txt"

# Directory to save the generated files
output_dir="${base_dir}/chloroplast/results/vg"
circ_output_dir="${output_dir}/circ"

# Ensure the output directories exist
mkdir -p "$output_dir"
mkdir -p "$circ_output_dir"

# Temporary file to store commands for parallel execution
commands_file=$(mktemp)

# Read each line from the file list and create commands
while IFS= read -r file; do
  # Remove the base directory from the path
  short_path="${file#$base_dir/}"

  # Extract the base name without the extension
  rank_name=$(basename "$short_path" .fasta.gz)

  # Define output file names
  vg_file="${output_dir}/${rank_name}.vg"
  stats_file="${output_dir}/${rank_name}.txt"
  circ_vg_file="${circ_output_dir}/${rank_name}_circ.vg"

  # Create commands to be run in parallel
  echo "vg construct -a -M <(zcat \"$file\") > \"$vg_file\" && vg stats -r \"$vg_file\" > \"$stats_file\" && vg circularize -a 2 -z \$(sed 's/[^,:]*://g' \"$stats_file\") \"$vg_file\" > \"$circ_vg_file\"" >> "$commands_file"
done < "$file_list"

# Run the commands in parallel
parallel --slf list_serv_2_14 < "$commands_file"

# Clean up the temporary commands file
rm "$commands_file"

echo "Individual file generation completed successfully."
