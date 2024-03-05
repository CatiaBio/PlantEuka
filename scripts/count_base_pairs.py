from Bio import SeqIO
import gzip
import sys

# Use the first command line argument as the gzipped fasta file path
fasta_file_path = sys.argv[1]

def count_bp_per_sequence(fasta_file):
    bp_counts = []
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            bp_counts.append((record.id, len(record.seq)))
    return bp_counts

# Call the function with the provided file path
bp_counts = count_bp_per_sequence(fasta_file_path)

# Print the counts
for id, count in bp_counts:
    print(f"{id}\t{count}")
