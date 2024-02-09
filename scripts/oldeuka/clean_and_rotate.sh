#!/bin/bash

mylist=last

cat $mylist | while read val;
do
    #FILE=/home/projects/animalia_mito/data/${val}/${val}_clean.fa.gz
    
    
    #if [[ -f $FILE ]]
    #then 
        #echo ${val}
        sed '/^>/!s/R/N/g; /^>/!s/Y/N/g; /^>/!s/W/N/g; /^>/!s/S/N/g; /^>/!s/K/N/g; /^>/!s/M/N/g; /^>/!s/D/N/g; /^>/!s/V/N/g; /^>/!s/H/N/g; /^>/!s/B/N/g; /^>/!s/X/N/g' <(zcat /home/projects/animalia_mito/data/$val/$val.fa.gz) > /home/projects/animalia_mito/data/${val}/${val}_clean.fa || sed '/^>/!s/R/N/g; /^>/!s/Y/N/g; /^>/!s/W/N/g; /^>/!s/S/N/g; /^>/!s/K/N/g; /^>/!s/M/N/g; /^>/!s/D/N/g; /^>/!s/V/N/g; /^>/!s/H/N/g; /^>/!s/B/N/g; /^>/!s/X/N/g' <(zcat /home/projects/animalia_mito/data/${val}/${val}_reps.fa.gz) > /home/projects/animalia_mito/data/${val}/${val}_clean.fa
        gzip -f /home/projects/animalia_mito/data/${val}/${val}_clean.fa
        nice -19 /home/ctools/opt/anaconda3_202105/bin/python3 /home/projects/animalia_mito/scripts/rotate_by_substring.py "${val}" > /home/projects/animalia_mito/data/${val}/rota_${val}.out
    #fi 
done

#echo "Process complete!" > /home/projects/animalia_mito/scripts/final.txt
