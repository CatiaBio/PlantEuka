import pandas as pd
import sys
import os
from collections import Counter

# General Description:
# This script processes a TSV file containing pairwise comparison data. 
# It filters pairs based on specified conditions, identifies unique accessions 
# to keep, removes duplicates, and counts the unique accessions before and after 
# filtering. The final results are saved to an output file.

# For Snakemake compatibility
file_path = snakemake.input.results_file
output_path = snakemake.output.filtered_file
organelle = snakemake.params.organelle

# Read the TSV file into a DataFrame
df = pd.read_csv(file_path, sep='\t')

# Remove duplicate pairs globally (considering both directions of pairs)
def normalize_pair(row):
    """Normalize pairs so (A,B) and (B,A) are treated as the same pair"""
    pair1, pair2 = row['pair1'], row['pair2']
    return tuple(sorted([pair1, pair2]))

# Add normalized pair column for deduplication
df['normalized_pair'] = df.apply(normalize_pair, axis=1)

# Remove duplicate pairs globally (keeping the first occurrence)
df = df.drop_duplicates(subset=['normalized_pair'], keep='first')

# Remove the temporary column
df = df.drop('normalized_pair', axis=1)

print(f"After deduplication: {len(df)} unique pairs remaining")

# Function to process pairs based on conditions
def process_pairs(df):
    """
    Process the input DataFrame to classify pairs into 'pairs_to_keep',
    'pairs_to_discard', and 'equal_pairs' based on specified conditions.
    
    Args:
    df (DataFrame): The input DataFrame containing pairwise comparison data.
    
    Returns:
    dict: A dictionary with keys as 'rank_name' and values as dictionaries
          containing lists of 'pairs_to_keep', 'pairs_to_discard', and 'equal_pairs'.
    """
    result = {}
    for rank_name, group in df.groupby(['rank', 'name']):
        rank, name = rank_name
        max_score = group['score'].max()
        threshold_score = max_score / 2
        
        # Initialize lists to store pairs
        pairs_to_keep = []
        pairs_to_discard = []
        equal_pairs = []
        
        # Filter pairs based on conditions
        for _, row in group.iterrows():
            pair = (row['pair1'], row['pair2'])
            if row['similarity'] > 99:
                equal_pairs.append(pair)
            elif (50 < row['similarity'] < 90) and (row['gaps'] < 15) and (row['score'] >= threshold_score):
                pairs_to_keep.append(pair)
            else:
                pairs_to_discard.append(pair)
        
        result[f'{rank}_{name}'] = {
            'pairs_to_keep': pairs_to_keep,
            'pairs_to_discard': pairs_to_discard,
            'equal_pairs': equal_pairs
        }
    
    return result

# Function to get unique accessions to keep
def get_unique_accessions(processed_pairs):
    """
    Determine unique accessions to keep for each 'rank_name' by removing duplicates
    based on 'equal_pairs' and count the unique accessions before and after filtering.
    
    Args:
    processed_pairs (dict): The dictionary returned from the process_pairs function.
    
    Returns:
    dict: A dictionary with 'rank_name' as keys and sets of unique accessions as values.
    dict: A dictionary with 'rank_name' as keys and tuples of (count_before, count_after) as values.
    list: A list of tuples containing original rank names and their corresponding counts.
    dict: A dictionary with 'rank_name' as keys and sets of discarded accessions as values.
    """
    final_result = {}
    counts = {}
    original_counts = []
    discarded_accessions = {}
    
    for rank_name, lists in processed_pairs.items():
        # Create a set of all accessions
        all_accessions = set()
        for pair in lists['pairs_to_keep']:
            all_accessions.update(pair)
        for pair in lists['pairs_to_discard']:
            all_accessions.update(pair)
        for pair in lists['equal_pairs']:
            all_accessions.update(pair)

        # Create a set of unique accessions from pairs_to_keep
        accessions_to_keep = set()
        for pair in lists['pairs_to_keep']:
            accessions_to_keep.update(pair)
        
        # Note: pairs_to_discard are already excluded from accessions_to_keep
        # as they are only collected from pairs_to_keep initially
        
        # Count unique accessions before filtering
        unique_accessions_before = set()
        for pair in lists['pairs_to_keep']:
            unique_accessions_before.update(pair)
        for pair in lists['pairs_to_discard']:
            unique_accessions_before.update(pair)
        for pair in lists['equal_pairs']:
            unique_accessions_before.update(pair)
        count_before = len(unique_accessions_before)
        
        # Sort equal_pairs by the first pair
        sorted_equal_pairs = sorted(lists['equal_pairs'], key=lambda x: x[0])
        
        # Iterate through equal_pairs to filter out duplicates
        # Remove the lexicographically larger accession ID to maintain consistency
        for pair in sorted_equal_pairs:
            pair1, pair2 = pair
            if pair1 in accessions_to_keep and pair2 in accessions_to_keep:
                # Remove the lexicographically larger accession to maintain consistency
                if pair1 > pair2:
                    accessions_to_keep.remove(pair1)
                else:
                    accessions_to_keep.remove(pair2)
        
        count_after = len(accessions_to_keep)
        counts[rank_name] = (count_before, count_after)
        final_result[rank_name] = accessions_to_keep
        
        # Determine discarded accessions
        discarded_accessions[rank_name] = all_accessions - accessions_to_keep
        
        original_counts.append((rank_name, count_before, count_after))
    
    return final_result, counts, original_counts, discarded_accessions

# Apply the function to process pairs and get unique accessions
processed_pairs = process_pairs(df)
unique_accessions, counts, original_counts, discarded_accessions = get_unique_accessions(processed_pairs)

# Save the result to a file 
os.makedirs(os.path.dirname(output_path), exist_ok=True)

with open(output_path, 'w') as f:
    for rank_name, accessions in unique_accessions.items():
        count_before, count_after = counts[rank_name]
        f.write(f'{rank_name} bf {count_before} af {count_after}\n')
        for accession in sorted(accessions):
            f.write(f'{accession}\n')
        f.write("Discarded\n")
        for discarded in sorted(discarded_accessions[rank_name]):
            f.write(f'{discarded}\n')

print(f'Results saved to {output_path}')
