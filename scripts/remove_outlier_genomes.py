import gzip

def remove_sequences(input_file, output_file, sequence_ids_to_remove):
    """
    Remove specific sequences from a gzipped FASTA file.

    Parameters:
    input_file (str): The file path to the input gzipped FASTA file.
    output_file (str): The file path to the output gzipped FASTA file.
    sequence_ids_to_remove (list): A list of IDs of sequences to be removed.
    """
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        write_sequence = True  # A flag to determine whether to write the sequence
        for line in infile:
            if line.startswith('>'):
                sequence_id = line.split()[0][1:]  # Assumes ID is the first part of the header line, without '>'
                # Check if the current sequence ID is in the list to remove
                if sequence_id in sequence_ids_to_remove:
                    write_sequence = False  # Do not write this sequence
                else:
                    write_sequence = True  # Write all other sequences
                    outfile.write(line)  # Write the header line
            else:
                if write_sequence:
                    outfile.write(line)  # Write the sequence line

# Example usage:
input_fasta = '../data/tmp_combined_sequences_mt.fa'  # Replace with your actual input file path
output_fasta = '../data/tmp_combined_sequences_clean_mt.fa'  # Replace with your desired output file path
sequence_ids = ['NC_031359.1', 'NC_031360.1']  # Replace with the IDs of the sequences you want to remove
remove_sequences(input_fasta, output_fasta, sequence_ids)
