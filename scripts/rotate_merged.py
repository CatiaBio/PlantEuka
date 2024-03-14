from Bio import SeqIO
import gzip
import os
import argparse

# Define the function to find the longest common substring
def findstem(arr):
    n = len(arr)
    s = arr[0]
    l = len(s)
    res = ""
    for i in range(l):
        for j in range(i + 1, l + 1):
            stem = s[i:j]
            if all(stem in arr[k] for k in range(1, n)):
                if len(res) < len(stem):
                    res = stem
    return res

# Define the function to rotate the sequences
def rotate(strg, n):
    return strg[n:] + strg[:n]

# Main function to process all FASTA files in the merged folder
def process_fasta_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith('.fa.gz'):
            file_path = os.path.join(directory, filename)
            seq_list = []
            
            # Read the sequences from the FASTA file
            with gzip.open(file_path, "rt") as file:
                for record in SeqIO.parse(file, "fasta"):
                    seq_list.append(str(record.seq))
            
            # Find the longest common substring
            common_substring = findstem(seq_list)
            
            # Rotate sequences if a common substring was found
            if common_substring:
                rotated_sequences = []
                with gzip.open(file_path, "rt") as file:
                    for record in SeqIO.parse(file, "fasta"):
                        index = str(record.seq).find(common_substring)
                        if index != -1:
                            record.seq = rotate(record.seq, index)
                        rotated_sequences.append(record)
                
                # Save the rotated sequences to a new file
                output_file_path = os.path.join(directory, f"{os.path.splitext(filename)[0]}_rotate.fa.gz")
                with gzip.open(output_file_path, "wt") as output_file:
                    SeqIO.write(rotated_sequences, output_file, "fasta")
                    
                print(f"Processed and saved: {output_file_path}")

# Set the directory path to the merged folder
merged_directory = "../data/chloroplast/merged"

# Call the main function to process the files
process_fasta_files(merged_directory)
