#!/bin/bash

for val in $(cat last);
do
    #FILE=/home/projects/animalia_mito/data/${val}/${val}_for_msa.fa.gz
    #if [ -f "$FILE" ]; then
    #    echo "Pairwise alignment for ${val} has already been performed. MSA file exists"
    #else
    nice -19 /home/ctools/opt/anaconda3_202105/bin/python3 /home/projects/animalia_mito/scripts/pairwise_needle.py ${val} > /home/projects/animalia_mito/data/${val}/pairwise_${val}.out
    nice -19 /home/ctools/opt/anaconda3_202105/bin/python3 /home/projects/animalia_mito/scripts/plot_matrix.py ${val} > /home/projects/animalia_mito/data/${val}/cleaned_${val}.out
    nice -19 zip -r -m /home/projects/animalia_mito/data/${val}/pw_scores.zip /home/projects/animalia_mito/data/${val}/*.txt -x /home/projects/animalia_mito/data/${val}/${val}_lin.txt
    
    #sleep 1m
    #fi

done
