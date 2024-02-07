import os
import argparse
import sys
import gzip
from Bio import SeqIO

def create_fam_fasta():
    os.system("ls /home/projects/animalia_mito/data/ > /home/projects/animalia_mito/scripts/all.txt")
    f = open("/home/projects/animalia_mito/scripts/all.txt", "r")
    f = f.read().splitlines()
    f = list(f)
    for i in f:
        if os.path.isfile(f"/home/projects/animalia_mito/data/{i}/{i}_reps.fa.gz"):
            pass
        else:
            os.system(f"zcat /home/projects/animalia_mito/data/{i}/*/*_reps.fa.gz /home/projects/animalia_mito/data/{i}/*/*/*_reps.fa.gz /home/projects/animalia_mito/data/{i}/*/*/*/*_reps.fa.gz /home/projects/animalia_mito/data/{i}/*/*/*/*/*_reps.fa.gz > /home/projects/animalia_mito/data/{i}/{i}.fa")
            os.system(f"gzip -f /home/projects/animalia_mito/data/{i}/{i}.fa")
            print(f"New multi-fasta file {i}.fa.gz was created.", file = sys.stdout)

def create_fam_lin():
    os.system("ls /home/projects/animalia_mito/data/ > /home/projects/animalia_mito/scripts/all.txt")
    f = open("/home/projects/animalia_mito/scripts/all.txt", "r")
    f = f.read().splitlines()
    f = list(f)
    for i in f:
        if os.path.isfile(f"/home/projects/animalia_mito/data/{i}/{i}_reps.fa.gz"):
            pass
        else:
            os.system(f"cat /home/projects/animalia_mito/data/{i}/*/*_lin.txt /home/projects/animalia_mito/data/{i}/*/*/*_lin.txt /home/projects/animalia_mito/data/{i}/*/*/*/*_lin.txt /home/projects/animalia_mito/data/{i}/*/*/*/*/*_lin.txt  > /home/projects/animalia_mito/data/{i}/{i}_lin.txt")
            print(f"New lineage file {i}_lin.txt was created.", file = sys.stdout)



new = create_fam_fasta()
lin = create_fam_lin()
