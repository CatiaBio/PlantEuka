#!/bin/bash

# Get the current working directory
base_dir=$(pwd)

# Directory where the generated files are saved
output_dir="${base_dir}/chloroplast/results/vg"
circ_output_dir="${output_dir}/circ"

# Ensure the output directories exist
mkdir -p "$output_dir"
mkdir -p "$circ_output_dir"

# Combine the generated rank_name_circ.vg files
circ_files=($(ls ${circ_output_dir}/*_circ.vg))
combined_circ_files=$(printf " %s" "${circ_files[@]}")
vg combine -c $combined_circ_files > "${output_dir}/combined_ranks_names.vg"
vg snarls "${output_dir}/combined_ranks_names.vg" > "${output_dir}/combined_ranks_names.txt"
vg view "${output_dir}/combined_ranks_names.vg" > "${output_dir}/combined_ranks_names.gfa"
vg gbwt -o "${output_dir}/combined_ranks_names.gbwt" -g "${output_dir}/combined_ranks_names.gg" -G "${output_dir}/combined_ranks_names.gfa"
vg index -j "${output_dir}/combined_ranks_names.dist" "${output_dir}/combined_ranks_names.vg"
vg minimizer -g "${output_dir}/combined_ranks_names.gbwt" -o "${output_dir}/combined_ranks_names.min" -k 20 -w 10 "${output_dir}/combined_ranks_names.vg"
vg convert -g "${output_dir}/combined_ranks_names.gfa" -o > "${output_dir}/combined_ranks_names.og"

echo "Combination of files and further processing completed successfully."
