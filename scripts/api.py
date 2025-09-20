import os
import joblib
import pandas as pd
from flask import Flask, request, jsonify
from google.cloud import storage
import io
import numpy as np

app = Flask(__name__)

# --- Global variables ---
model = None
storage_client = None
bucket = None
feature_columns = None

def download_model_and_features():
    """Downloads the model and feature columns from GCS."""
    global model, feature_columns
    print("Downloading model and features from GCS...")

    # Download model
    model_blob = bucket.blob("model.joblib")
    model_data = model_blob.download_as_string()
    model = joblib.load(io.BytesIO(model_data))
    print("Model downloaded and loaded successfully.")

    # Download feature columns
    features_blob = bucket.blob("features.csv")
    features_data = features_blob.download_as_string()
    features_df = pd.read_csv(io.StringIO(features_data.decode("utf-8")), nrows=0) # Read only headers
    feature_columns = features_df.columns.tolist()
    print(f"Feature columns loaded. Expecting {len(feature_columns)} features.")


@app.before_request
def load_model_and_features():
    """Load the model and feature columns before the first request."""
    global storage_client, bucket
    if model is None:
        storage_client = storage.Client()
        bucket = storage_client.get_bucket("zntic-data")
        download_model_and_features()

@app.route('/predict', methods=['POST'])
def predict():
    """Receives a request and returns a prediction."""
    if not request.json or 'features' not in request.json:
        return jsonify({'error': 'Invalid input: JSON with "features" key is required.'}), 400

    try:
        features = request.json['features']
        
        if len(features) != len(feature_columns):
            return jsonify({'error': f'Invalid number of features. Expected {len(feature_columns)}, got {len(features)}.'}), 400

        # Create a DataFrame with the correct column names and order
        df = pd.DataFrame([features], columns=feature_columns)
        df = df[feature_columns] # Ensure column order

        # Make a prediction
        prediction = model.predict(df.drop('zoonotic', axis=1))
        prediction_proba = model.predict_proba(df.drop('zoonotic', axis=1))

        # Return the result
        result = {
            'prediction': 'zoonotic' if prediction[0] == 1 else 'non-zoonotic',
            'confidence': {
                'non-zoonotic': prediction_proba[0][0],
                'zoonotic': prediction_proba[0][1]
            }
        }
        return jsonify(result)

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return "OK", 200

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))