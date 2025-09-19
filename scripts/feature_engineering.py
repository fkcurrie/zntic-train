import os
import itertools
from Bio import SeqIO
import pandas as pd

# Input and output file paths
input_file = os.path.join("data", "processed_genomes.fna")
output_file = os.path.join("data", "kmer_features.csv")

# K-mer size
k = 4

# Check if the input data file exists
if not os.path.exists(input_file):
    print(f"Error: {input_file} not found. Please run the preprocessing script first.")
    exit()

# Generate all possible k-mers
bases = ['A', 'C', 'G', 'T']
all_kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
kmer_dict = {kmer: 0 for kmer in all_kmers}

# List to store feature dictionaries for each sequence
features_list = []

# Process sequences in batches to manage memory
batch_size = 1000
record_iterator = SeqIO.parse(input_file, "fasta")

while True:
    batch = list(itertools.islice(record_iterator, batch_size))
    if not batch:
        break

    for record in batch:
        # Initialize k-mer counts for the current sequence
        current_kmer_counts = kmer_dict.copy()
        
        # Count k-mers in the sequence
        sequence = str(record.seq).upper()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer in current_kmer_counts:
                current_kmer_counts[kmer] += 1
        
        # Normalize the counts to get frequencies
        total_kmers = len(sequence) - k + 1
        if total_kmers > 0:
            for kmer in current_kmer_counts:
                current_kmer_counts[kmer] /= total_kmers
        
        # Add sequence ID and features to the list
        current_kmer_counts['sequence_id'] = record.id
        features_list.append(current_kmer_counts)

    print(f"Processed {len(features_list)} sequences...")

# Create a DataFrame and save to CSV
df = pd.DataFrame(features_list)
# Reorder columns to have sequence_id first
cols = ['sequence_id'] + [col for col in df.columns if col != 'sequence_id']
df = df[cols]
df.to_csv(output_file, index=False)

print(f"Feature engineering complete. Features saved to {output_file}")
