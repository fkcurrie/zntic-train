import os
from google.cloud import storage
from Bio import Entrez

# Configuration
BUCKET_NAME = "zntic-train-data"
DESTINATION_BLOB_NAME = "raw_data/avian_influenza.fna"
Entrez.email = "frank@sfle.ca"

def download_and_upload_to_gcs():
    """Downloads influenza data from NCBI and uploads it to GCS."""

    print("Searching NCBI for avian influenza genomes...")
    # Search for avian influenza complete genomes
    handle = Entrez.esearch(
        db="nucleotide",
        term="Influenza A virus (H5N1) complete genome",
        retmax=500  # Using a smaller number for a quicker initial run
    )
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    print(f"Found {len(id_list)} records.")

    if not id_list:
        print("No records found. Exiting.")
        return

    print("Fetching records from NCBI...")
    # Fetch the records in FASTA format
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
    fasta_records = handle.read()
    handle.close()

    print(f"Uploading {len(fasta_records)} bytes to GCS bucket '{BUCKET_NAME}'...")
    # Upload to GCS
    storage_client = storage.Client()
    bucket = storage_client.bucket(BUCKET_NAME)
    blob = bucket.blob(DESTINATION_BLOB_NAME)

    blob.upload_from_string(fasta_records)

    print(f"File successfully uploaded to gs://{BUCKET_NAME}/{DESTINATION_BLOB_NAME}")

if __name__ == "__main__":
    download_and_upload_to_gcs()
