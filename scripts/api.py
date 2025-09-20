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

def download_model():
    """Downloads the model from GCS."""
    global model
    print("Downloading model from GCS...")
    blob = bucket.blob("model.joblib")
    model_data = blob.download_as_string()
    model = joblib.load(io.BytesIO(model_data))
    print("Model downloaded and loaded successfully.")

@app.before_request
def load_model():
    """Load the model before the first request."""
    global storage_client, bucket
    if model is None:
        storage_client = storage.Client()
        bucket = storage_client.get_bucket("zntic-data")
        download_model()

@app.route('/predict', methods=['POST'])
def predict():
    """Receives a request and returns a prediction."""
    if not request.json or 'features' not in request.json:
        return jsonify({'error': 'Invalid input: JSON with "features" key is required.'}), 400

    try:
        # Create a DataFrame from the input features
        # The model expects a DataFrame with specific column names.
        # For this placeholder, we assume the input is a list of feature values
        # and we'll create a DataFrame with dummy column names.
        features = request.json['features']
        
        # Get the expected number of features from the model
        if hasattr(model, 'n_features_in_'):
            num_features = model.n_features_in_
        else:
            # Fallback for older scikit-learn versions or different models
            return jsonify({'error': 'Could not determine the number of features expected by the model.'}), 500

        if len(features) != num_features:
            return jsonify({'error': f'Invalid number of features. Expected {num_features}, got {len(features)}.'}), 400

        # Create a DataFrame with the correct number of columns
        df = pd.DataFrame([features], columns=[f'feature_{i}' for i in range(num_features)])

        # Make a prediction
        prediction = model.predict(df)
        prediction_proba = model.predict_proba(df)

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
