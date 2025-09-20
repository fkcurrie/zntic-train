import os
import pandas as pd
from google.cloud import storage
import io
import joblib
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import numpy as np

def main():
    """
    Main function to train the model.
    """
    # Set up GCS client
    storage_client = storage.Client()
    bucket = storage_client.get_bucket("zntic-data")

    # Download the features from GCS
    print("Downloading features from GCS...")
    blob = bucket.blob("features.csv")
    features_data = blob.download_as_string()

    # Load the features into a pandas DataFrame
    print("Loading features into DataFrame...")
    df = pd.read_csv(io.StringIO(features_data.decode("utf-8")))

    # --- Placeholder for Labels ---
    # In a real scenario, you would have a separate file with labels.
    # For now, we'll generate a placeholder 'zoonotic' column.
    print("Generating placeholder labels...")
    np.random.seed(42) # for reproducibility
    df['zoonotic'] = np.random.randint(0, 2, df.shape[0])
    # --- End Placeholder ---

    # Split the data into training and testing sets
    print("Splitting data into training and testing sets...")
    X = df.drop('zoonotic', axis=1)
    y = df['zoonotic']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train a simple Logistic Regression model
    print("Training Logistic Regression model...")
    model = LogisticRegression(max_iter=1000)
    model.fit(X_train, y_train)

    # Evaluate the model
    print("Evaluating model...")
    y_pred = model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    print(f"Model Accuracy: {acc:.2f}")

    # Save the trained model to a file
    print("Saving model to file...")
    model_filename = 'model.joblib'
    joblib.dump(model, model_filename)

    # Upload the model to GCS
    print("Uploading model to GCS...")
    model_blob = bucket.blob(model_filename)
    model_blob.upload_from_filename(model_filename)

    print("Training complete. Model uploaded to GCS.")

if __name__ == "__main__":
    main()
