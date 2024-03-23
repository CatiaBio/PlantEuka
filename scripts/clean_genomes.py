import gzip
import os
from datetime import datetime
import sys
from collections import Counter

def clean_fasta(fasta_path):
    with gzip.open(fasta_path, 'rt') as file:
        original_content = file.read()
    replaced_characters = Counter()
    cleaned_lines = []
    for line in original_content.split('\n'):
        if line.startswith('>'):
            cleaned_lines.append(line)
        else:
            cleaned_line, counts = clean_sequence(line)
            cleaned_lines.append(cleaned_line)
            replaced_characters.update(counts)
    cleaned_content = '\n'.join(cleaned_lines)
    return cleaned_content, replaced_characters

def clean_sequence(sequence):
    replaced_characters = Counter()
    cleaned_sequence = []
    for nucleotide in sequence:
        if nucleotide not in 'ACGTN':
            cleaned_sequence.append('N')
            replaced_characters[nucleotide] += 1
        else:
            cleaned_sequence.append(nucleotide)
    return ''.join(cleaned_sequence), replaced_characters

def process_directory(base_dir, log_file_name):
    categories = ["genus", "family", "order"]
    start_time = datetime.now()

    with open(log_file_name, 'w') as log_file:
        log_file.write(f"Starting cleanup process at {start_time}\n")

        for category in categories:
            category_path = os.path.join(base_dir, category)
            original_before_clean_path = os.path.join(category_path, 'original_before_clean')
            changes_log_path = os.path.join(original_before_clean_path, 'changes_log.tsv')

            if not os.path.isdir(category_path):
                log_file.write(f"Directory does not exist: {category_path}\n")
                continue

            os.makedirs(original_before_clean_path, exist_ok=True)

            with open(changes_log_path, 'w') as changes_log:
                changes_log.write("ID\tNumber_Replacements\tReplaced_Characters\n")

                for root, dirs, files in os.walk(category_path):
                    for file in files:
                        if file.endswith('.fasta.gz'):
                            fasta_path = os.path.join(root, file)
                            base_name = os.path.splitext(os.path.splitext(os.path.basename(fasta_path))[0])[0]

                            cleaned_content, replaced_characters = clean_fasta(fasta_path)
                            if replaced_characters:
                                cleaned_path = fasta_path.replace('.fasta.gz', '_cleaned.fasta.gz')
                                with gzip.open(cleaned_path, 'wt') as cleaned_file:
                                    cleaned_file.write(cleaned_content)

                                original_backup_path = os.path.join(original_before_clean_path, file)
                                os.rename(fasta_path, original_backup_path)
                                os.rename(cleaned_path, fasta_path)

                                replaced_info = ', '.join(f"{char}({count})" for char, count in sorted(replaced_characters.items()))
                                changes_log.write(f"{base_name}\t{sum(replaced_characters.values())}\t{replaced_info}\n")
                                log_file.write(f"Cleaned and replaced FASTA file: {fasta_path}\n")

        log_file.write(f"All analyses completed at {datetime.now()}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <base_directory> <log_file_name>")
        sys.exit(1)

    base_dir = sys.argv[1]
    log_file_name = sys.argv[2]

    process_directory(base_dir, log_file_name)
