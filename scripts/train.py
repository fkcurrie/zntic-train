import os
import pandas as pd
from google.cloud import storage
import io

def main():
    """
    Main function to train the model.
    """
    # Set up GCS client
    storage_client = storage.Client()
    bucket = storage_client.get_bucket("zntic-data")

    # Download the features from GCS
    blob = bucket.blob("features.csv")
    features_data = blob.download_as_string()

    # Load the features into a pandas DataFrame
    df = pd.read_csv(io.StringIO(features_data.decode("utf-8")))

    # Print the head of the DataFrame
    print(df.head())

if __name__ == "__main__":
    main()