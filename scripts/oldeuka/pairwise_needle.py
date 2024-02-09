from itertools import combinations 
from Bio import SeqIO
import gzip 
import os
import pandas as pd
import re
from time import sleep
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('path', help= "name list")
args = parser.parse_args()

need_needle = []
try:
    path = f"/home/projects/animalia_mito/data/{args.path}/{args.path}_rota_new.fa.gz"

    with gzip.open(path, "rt") as file:
        multi = combinations(SeqIO.parse(file, "fasta"), 2)
        for pair in multi:
            seq1 = pair[0].seq
            seq2 = pair[1].seq
            if os.path.isfile(f"/home/projects/animalia_mito/data/{args.path}/{args.path}_{pair[0].id}_{pair[1].id}.txt") and os.path.getsize(f"/home/projects/animalia_mito/data/{args.path}/{args.path}_{pair[0].id}_{pair[1].id}.txt") > 0:
                pass
            else:
                os.system(f"/home/ctools/EMBOSS-6.6.0/emboss/needle -asequence asis:{seq1} -bsequence asis:{seq2} -outfile /home/projects/animalia_mito/data/{args.path}/{args.path}_{pair[0].id}_{pair[1].id}.txt -gapopen 10 -gapextend 3")
                sleep(10)
except:
    print(f"{args.path}", file = sys.stdout)


        























        #f = open (f"../data/{args.path}/{args.path}_{pair[0].id}_{pair[1].id}.txt", "r")
        #f = f.read().split('\n')
        #for line in f:
           #if re.search("Score", line):
               #l = float(line[9:])
               #df.at[pair[1].id,pair[0].id] = l
               #df.at[pair[0].id,pair[1].id] = l

#print(df)
#x = []
#for i in name:
    #print(i)
    #x1 = df[i].sum()
    #print(x1)
    #x.append(x1)

#print(x)

#plt.figure(figsize=(20,10))
#plt.bar(name, x)
#plt.xticks(rotation='vertical')
#plt.title(f"Pairwise alignment score sums for {args.path}")
#plt.savefig(f"../data/{args.path}/{args.path}_pairwise_sum.png")
