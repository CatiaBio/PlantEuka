from itertools import combinations
from Bio import SeqIO
import gzip
import pandas as pd
import sys
import re
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('path', help= "Path to file list")
args = parser.parse_args()

path = f"/home/projects/MAAG/PlantEuka/genomes/chloroplast/original/*.fasta.gz"
#path = f"/home/projects/animalia_mito/data/{args.path}/{args.path}_rota_new.fa.gz"

# Load the data from the TSV file
file_path = 'pairwise_results_cp.tsv'
df = pd.read_csv(file_path, sep='\t')
output_file = "genomes_to_exclude.txt"


# Group the data by 'Name'
grouped = df.groupby('Name')

score = []
similarity = []
gaps = []
duplicate_samples = []

for name, group in grouped:
    max_score = group['Score'].max()
    half_max_score = max_score / 2
    
    for entry_index, entry in group.iterrows():

        if (entry['Similarity'] < 49.9):
            similarity.append(entry['pair1'],entry['pair2'])
        if (entry['Gaps'] > 15):
            gaps.append(entry['pair1'],entry['pair2'])
        if (entry['Score'] < half_max_score):
            score.append(entry['pair1'],entry['pair2'])
        if (entry['Similarity'] > 99):
            duplicate_samples.append(entry['pair1'])


all_exclude = []

if len(duplicate_samples) != 0:
    all_exclude.append(duplicate_samples)
    print("The following sample will be excluded as it is a duplicate sample")
else:
    pass
if len(score) != 0:
    all_exclude.append(score)
    print("The following samples will be excluded due to a low pairwise score")
else:
    pass
if len(similarity) != 0:
    all_exclude.append(similarity)
    print("The following samples will be excluded due to a low similarity score")
else:
    pass
if len(gaps) != 0:
    all_exclude.append(gaps)
    print("The following samples will be excluded due to a high gap score")
else:
    pass

for item in all_exclude:
    print(item)


