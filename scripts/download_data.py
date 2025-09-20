import os
from Bio import Entrez
from Bio import SeqIO
import subprocess

def download_and_upload_to_gcs():
    """
    Downloads avian influenza data from NCBI and uploads it to a GCS bucket.
    """
    # Set your email for NCBI Entrez
    Entrez.email = "user@example.com"

    # Search for avian influenza complete genomes
    print("Searching NCBI for avian influenza genomes...")
    handle = Entrez.esearch(db="nucleotide", term="Influenza A virus (H5N1) complete genome", retmax=500)
    record = Entrez.read(handle)
    handle.close()
    print(f"Found {len(record['IdList'])} records.")

    # Fetch the records in FASTA format
    print("Fetching records from NCBI...")
    id_list = record["IdList"]
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
    fasta_records = handle.read()
    handle.close()

    # Save the records to a local file
    local_file = "avian_influenza.fasta"
    with open(local_file, "w") as f:
        f.write(fasta_records)

    # Upload to GCS
    print(f"Uploading {local_file} to GCS bucket 'zntic-data'...")
    subprocess.run(["gcloud", "storage", "cp", local_file, "gs://zntic-data/"], check=True)

    print("Upload complete.")

if __name__ == "__main__":
    download_and_upload_to_gcs()