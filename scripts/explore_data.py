import os
from Bio import SeqIO

# Path to the downloaded data
data_file = os.path.join("data", "ncbi_dataset", "data", "genomic.fna")

# Check if the data file exists
if not os.path.exists(data_file):
    print(f"Error: {data_file} not found. Please run download_data.py first.")
    exit()

# Parse the FASTA file
records = list(SeqIO.parse(data_file, "fasta"))

# Print the number of records
print(f"Number of records: {len(records)}")

# Print the ID, description, and length of each record
for record in records:
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Length: {len(record.seq)}")
    print("-" * 20)
