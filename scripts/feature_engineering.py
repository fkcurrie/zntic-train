import os
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from google.cloud import storage
import io

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
    # Set up GCS client
    storage_client = storage.Client()
    bucket = storage_client.get_bucket("zntic-data")

    # Download the data from GCS
    blob = bucket.blob("avian_influenza.fasta")
    fasta_records = blob.download_as_string().decode("utf-8")

    # Parse the FASTA file
    records = list(SeqIO.parse(io.StringIO(fasta_records), "fasta"))

    # Calculate k-mer composition for each record
    k = 3
    kmer_compositions = []
    for record in records:
        kmer_composition = get_kmer_composition(str(record.seq), k)
        kmer_compositions.append(kmer_composition)

    # Create a pandas DataFrame from the k-mer compositions
    df = pd.DataFrame(kmer_compositions)
    df = df.fillna(0)

    # Save the DataFrame to a CSV file in GCS
    output_blob = bucket.blob("features.csv")
    output_blob.upload_from_string(df.to_csv(index=False), "text/csv")

    print("Saved features to GCS")

if __name__ == "__main__":
    main()
