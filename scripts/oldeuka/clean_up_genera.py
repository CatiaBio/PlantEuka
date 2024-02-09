import os
import gzip
import argparse
from Bio import SeqIO

p = argparse.ArgumentParser()
p.add_argument('path')
args = p.parse_args()


artho = open(args.path, "r")
arth_new = artho.read().splitlines()
artho_new = list(arth_new)


for fam in artho_new:
    if os.path.isfile(f"../data/{fam}/{fam}_reps.fa.gz"):
        new_tx = []
        exclude = []
        with open(f"../data/{fam}/{fam}_lin.txt") as lin:

            lines = lin.readlines()
            for line in lines:

                split = line.split(':')
                accession_number = split[0]
                list_names = [x.strip().replace('\'','') for x in split[1].lstrip('[').rstrip(']\n').split(',')]
                new_fam = fam.lower()

                if any([x.endswith(f'{new_fam}') for x in list_names]):
                    exclude.append(accession_number)
                else:
                    new_tx.append(line)
                    

                
                
        
        with gzip.open(f"/home/projects/animalia_mito/data/{fam}/{fam}_reps.fa.gz", "rt") as fam_file:
            seq = SeqIO.parse(fam_file, "fasta")
            seqs = []
            for m in seq:
                if m.id not in exclude:
                    seqs.append(m)
                else:
                    pass
        with gzip.open(f"../data/{fam}/{fam}_reps.fa.gz", "wt") as g:
            SeqIO.write(seqs, g, "fasta")

        with open(f"../data/{fam}/{fam}_lin.txt", "wt") as k:
            for item in new_tx:
                k.write(item)


    
        

