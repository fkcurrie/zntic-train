import os
from Bio import SeqIO

# Input and output file paths
input_file = os.path.join("data", "ncbi_dataset", "data", "genomic.fna")
output_file = os.path.join("data", "processed_genomes.fna")

# Minimum sequence length threshold
min_length = 1000

# Check if the input data file exists
if not os.path.exists(input_file):
    print(f"Error: {input_file} not found. Please run the download script first.")
    exit()

# Filter sequences and write to the output file
with open(output_file, "w") as out_handle:
    count = 0
    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq) >= min_length:
            SeqIO.write(record, out_handle, "fasta")
            count += 1

print(f"Processed {count} records with length >= {min_length} bp.")
print(f"Filtered genomes saved to {output_file}")
