#!/bin/bash

for val in $(cat euka_db); 
do
    #FILE=/home/projects/animalia_mito/data/${val}/circ_${val}.vg

    #if [ -s "$FILE" ]; then
    #    echo "The circ_${val}.vg file exists and is not empty. Graph construction for ${val} has already been performed."
    #else

     #   nice -19 /home/ctools/vg/bin/vg construct -a -M /home/projects/animalia_mito/data/${val}/${val}_prank.best.fas > /home/projects/animalia_mito/data/${val}/${val}.vg
        nice -19 /home/ctools/vg/bin/vg stats -r /home/projects/animalia_mito/data/${val}/${val}.vg > /home/projects/animalia_mito/data/${val}/pos.txt
        end_node=$(sed 's/[^,:]*://g' /home/projects/animalia_mito/data/${val}/pos.txt)
        nice -19 /home/ctools/vg/bin/vg circularize -a 2 -z $end_node /home/projects/animalia_mito/data/${val}/${val}.vg > /home/projects/animalia_mito/data/${val}/circ_${val}.vg
    #fi
done

echo "Process complete!" > /home/projects/animalia_mito/scripts/final.txt
