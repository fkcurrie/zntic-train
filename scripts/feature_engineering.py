import os
import itertools
from Bio import SeqIO
import pandas as pd
from google.cloud import storage

# Configuration
BUCKET_NAME = "zntic-train-data"
SOURCE_BLOB_NAME = "raw_data/avian_influenza.fna"
DESTINATION_BLOB_NAME = "features/kmer_features.csv"
K = 4

def generate_kmer_features():
    """Downloads genome data from GCS, calculates k-mer frequencies, and uploads the results."""

    storage_client = storage.Client()
    bucket = storage_client.bucket(BUCKET_NAME)

    # Download the data from GCS
    print(f"Downloading gs://{BUCKET_NAME}/{SOURCE_BLOB_NAME}...")
    blob = bucket.blob(SOURCE_BLOB_NAME)
    fasta_string = blob.download_as_string().decode('utf-8')
    
    # Generate all possible k-mers
    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in itertools.product(bases, repeat=K)]
    kmer_dict = {kmer: 0 for kmer in all_kmers}

    features_list = []
    
    # Use a file-like object to parse the FASTA string
    from io import StringIO
    fasta_io = StringIO(fasta_string)

    print("Calculating k-mer frequencies...")
    for record in SeqIO.parse(fasta_io, "fasta"):
        current_kmer_counts = kmer_dict.copy()
        sequence = str(record.seq).upper()
        
        # Count k-mers
        for i in range(len(sequence) - K + 1):
            kmer = sequence[i:i+K]
            if kmer in current_kmer_counts:
                current_kmer_counts[kmer] += 1
        
        # Normalize to get frequencies
        total_kmers = len(sequence) - K + 1
        if total_kmers > 0:
            for kmer in current_kmer_counts:
                current_kmer_counts[kmer] /= total_kmers
        
        current_kmer_counts['sequence_id'] = record.id
        features_list.append(current_kmer_counts)

    print(f"Processed {len(features_list)} sequences.")

    # Create a DataFrame and save to CSV in memory
    df = pd.DataFrame(features_list)
    cols = ['sequence_id'] + [col for col in df.columns if col != 'sequence_id']
    df = df[cols]
    
    output_csv = df.to_csv(index=False)

    # Upload the features to GCS
    print(f"Uploading features to gs://{BUCKET_NAME}/{DESTINATION_BLOB_NAME}...")
    feature_blob = bucket.blob(DESTINATION_BLOB_NAME)
    feature_blob.upload_from_string(output_csv)
    
    print("Feature engineering complete.")

if __name__ == "__main__":
    generate_kmer_features()
