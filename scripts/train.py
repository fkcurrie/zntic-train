import os
import pandas as pd

def main():
    """
    Main function to train the model.
    """
    # Path to the features file
    features_file = os.path.join("data", "features.csv")

    # Check if the features file exists
    if not os.path.exists(features_file):
        print(f"Error: {features_file} not found. Please run feature_engineering.py first.")
        exit()

    # Load the features
    df = pd.read_csv(features_file)

    # Print the head of the DataFrame
    print(df.head())

if __name__ == "__main__":
    main()
