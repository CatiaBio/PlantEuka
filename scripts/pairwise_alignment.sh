#!/bin/bash

# Root directory for chloroplast sorted sequences
root_dir="genomes/chloroplast/sorted"

# Output directory for pairwise alignments, change this to a central location if needed
output_root_dir="genomes/chloroplast/pairwise_distances"
mkdir -p "$output_root_dir"

gapopen=16
gapextend=4

# Iterate over genus, family, and order directories
for taxonomy in genus family order; do
    # Iterate over subdirectories within each taxonomy
    for species_dir in "$root_dir"/$taxonomy/*; do
        if [ -d "$species_dir" ]; then
            species=$(basename "$species_dir")
            output_dir="${output_root_dir}/${taxonomy}_${species}"
            mkdir -p "$output_dir"
            
            # Decompress .fasta.gz files for processing
            find "$species_dir" -type f -name "*.fasta.gz" | while read fasta_gz; do
                fasta="${fasta_gz%.fasta.gz}.fa"
                gunzip -c "$fasta_gz" > "$fasta"
                echo "Decompressed $fasta"
            done
            
            # Perform pairwise alignments
            find "$species_dir" -type f -name "*.fa" | while read seq1; do
                find "$species_dir" -type f -name "*.fa" | while read seq2; do
                    if [ "$seq1" != "$seq2" ]; then
                        base1=$(basename "$seq1" .fa)
                        base2=$(basename "$seq2" .fa)
                        stretcher_outfile="${output_dir}/${base1}_vs_${base2}.needle"
                        echo "Aligning $base1 vs $base2"
                        /home/ctools/EMBOSS-6.6.0/emboss/stretcher -asequence "$seq1" -bsequence "$seq2" -gapopen $gapopen -gapextend $gapextend -outfile "$stretcher_outfile"
                    fi
                done
            done
            
            # Cleanup decompressed .fa files
            find "$species_dir" -type f -name "*.fa" -exec rm {} +
        fi
    done
done

echo "Pairwise alignment process completed."

