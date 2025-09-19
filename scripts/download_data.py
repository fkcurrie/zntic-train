import os
from Bio import Entrez
from Bio import SeqIO

# Set your email for NCBI Entrez
Entrez.email = "user@example.com"

# Create the data directory if it doesn't exist
if not os.path.exists("data"):
    os.makedirs("data")

# Search for avian influenza complete genomes
handle = Entrez.esearch(db="nucleotide", term="Influenza A virus (H5N1) complete genome", retmax=100)
record = Entrez.read(handle)
handle.close()

# Fetch the records in FASTA format
id_list = record["IdList"]
handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
fasta_records = handle.read()
handle.close()

# Save the records to a file
output_file = os.path.join("data", "avian_influenza.fasta")
with open(output_file, "w") as f:
    f.write(fasta_records)

print(f"Downloaded {len(id_list)} records to {output_file}")
