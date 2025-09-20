import os
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

def get_kmer_composition(sequence, k):
    """
    Calculates the k-mer composition of a DNA sequence.
    """
    kmer_counts = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] += 1
    return kmer_counts

def main():
    """
    Main function to perform feature engineering.
    """
    # Path to the downloaded data
    data_file = os.path.join("data", "ncbi_dataset", "data", "genomic.fna")

    # Check if the data file exists
    if not os.path.exists(data_file):
        print(f"Error: {data_file} not found. Please run download_data.py first.")
        exit()

    # Parse the FASTA file
    records = list(SeqIO.parse(data_file, "fasta"))

    # Calculate k-mer composition for each record
    k = 3
    kmer_compositions = []
    for record in records:
        kmer_composition = get_kmer_composition(str(record.seq), k)
        kmer_compositions.append(kmer_composition)

    # Create a pandas DataFrame from the k-mer compositions
    df = pd.DataFrame(kmer_compositions)
    df = df.fillna(0)

    # Save the DataFrame to a CSV file
    output_file = os.path.join("data", "features.csv")
    df.to_csv(output_file, index=False)

    print(f"Saved features to {output_file}")

if __name__ == "__main__":
    main()