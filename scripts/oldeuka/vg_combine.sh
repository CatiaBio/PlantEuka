#!/bin/bash

touch all_clades
for item in $(cat euka_db);
do
    echo "/home/projects/animalia_mito/data/${item}/circ_${item}.vg" >> all_clades
done

tr '\n' ' ' < all_clades > new_all_clades

