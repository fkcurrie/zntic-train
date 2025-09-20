import os
from Bio import Entrez
from Bio import SeqIO
from google.cloud import storage
from google.oauth2 import service_account

def download_and_upload_to_gcs():
    """
    Downloads avian influenza data from NCBI and uploads it to a GCS bucket.
    """
    # Set up GCS client with service account
    credentials = service_account.Credentials.from_service_account_file(
        "/app/secrets/gcr-key.json"
    )
    storage_client = storage.Client(credentials=credentials)
    bucket = storage_client.get_bucket("zntic-data")

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

    # Upload to GCS
    print(f"Uploading {len(fasta_records)} bytes to GCS bucket 'zntic-data'...")
    blob = bucket.blob("avian_influenza.fasta")
    blob.upload_from_string(fasta_records)

    print("Upload complete.")

if __name__ == "__main__":
    download_and_upload_to_gcs()
