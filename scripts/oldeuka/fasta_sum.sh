#!/bin/bash

mylist=really_all.txt
touch db_seqs.fa
cat $mylist | while read val;
do
    zcat ../data/${val}/${val}_reps.fa.gz >> db_seqs.fa
done

gzip db_seqs.fa 
