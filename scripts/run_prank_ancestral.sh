#!/bin/bash

for val in $(cat xaq);
do
    file=/home/projects/animalia_mito/data/${val}/ancestral_reconstruction/${val}_prank.best.fas
    if [ -s "$file" ]; then
        echo "Prank MSA has been run already for ${val} clade. Proceeding to the next clade."
    else
        zcat /home/projects/animalia_mito/data/${val}/${val}_for_msa.fa.gz | nice -19 /home/ctools/prank-msa/src/prank -d=/dev/stdin -o=/home/projects/animalia_mito/data/${val}/${val}_prank -DNA
    fi
done
