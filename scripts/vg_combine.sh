#!/bin/bash

# Get the current working directory
base_dir=$(pwd)

# Directory where the generated files are saved
circularize_dir="${base_dir}/chloroplast/results/vg/circularize"
combine_output_dir="${base_dir}/chloroplast/results/vg/combine"

# Ensure the output directory exists
mkdir -p "$combine_output_dir"

# Combine the generated rank_name_circ.vg files
circ_files=($(ls ${circularize_dir}/*_circ.vg))
combined_circ_files=$(printf " %s" "${circ_files[@]}")
vg combine -c $combined_circ_files > "${combine_output_dir}/planteuka_db.vg"
vg snarls "${combine_output_dir}/planteuka_db.vg" > "${combine_output_dir}/planteuka_db.txt"
vg view "${combine_output_dir}/planteuka_db.vg" > "${combine_output_dir}/planteuka_db.gfa"
vg gbwt -o "${combine_output_dir}/planteuka_db.gbwt" -g "${combine_output_dir}/planteuka_db.gg" -G "${combine_output_dir}/planteuka_db.gfa"
vg index -j "${combine_output_dir}/planteuka_db.dist" "${combine_output_dir}/planteuka_db.vg"
vg minimizer -g "${combine_output_dir}/planteuka_db.gbwt" -o "${combine_output_dir}/planteuka_db.min" -k 20 -w 10 "${combine_output_dir}/planteuka_db.vg"
vg convert -g "${combine_output_dir}/planteuka_db.gfa" -o > "${combine_output_dir}/planteuka_db.og"

echo "Combination of files and further processing completed successfully."