from Bio import SeqIO
import gzip
import argparse

p = argparse.ArgumentParser()
p.add_argument('family')
args = p.parse_args()

#https://www.geeksforgeeks.org/longest-common-substring-array-strings/
def findstem(arr):
 
    # Determine size of the array
    n = len(arr)
 
    # Take first word from array
    # as reference
    s = arr[0]
    l = len(s)
 
    res = ""
 
    for i in range(l):
        for j in range(i + 1, l + 1):
 
            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = s[i:j]
            k = 1
            for k in range(1, n):
 
                # Check if the generated stem is
                # common to all words
                if stem not in arr[k]:
                    break
 
            # If current substring is present in
            # all strings and its length is greater
            # than current result
            if (k + 1 == n and len(res) < len(stem)):
                res = stem
 
    return res


def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr

seq_list = []

with gzip.open(f"/home/projects/animalia_mito/data/{args.family}/{args.family}_clean.fa.gz", "rt") as file:
#with gzip.open("/home/projects/animalia_mito/scripts/test.fa.gz", "rt") as file:
    multi = SeqIO.parse(file, "fasta")
    for record in multi:
        seq_list.append(str(record.seq))

    
print("The sequence list has length: ")
print(len(seq_list))
#x = long_substr(seq_list)
x = findstem(seq_list)
print(x)

def rotate(strg, n):
    # Rotate the string by n
    return strg[n:] + strg[:n]

rota = []
with gzip.open(f"/home/projects/animalia_mito/data/{args.family}/{args.family}_clean.fa.gz", "rt") as f:
#with gzip.open("/home/projects/animalia_mito/scripts/test.fa.gz", "rt") as f:
    msa = SeqIO.parse(f, "fasta")
    
    for record in msa:
        record.seq = rotate(record.seq, record.seq.find(str(x)))
        rota.append(record)
        
    with gzip.open(f"/home/projects/animalia_mito/data/{args.family}/{args.family}_rota_new.fa.gz", "wt") as g:
    #with gzip.open("/home/projects/animalia_mito/scripts/test_rota.fa.gz", "wt") as g:
        SeqIO.write(rota, g, "fasta")
