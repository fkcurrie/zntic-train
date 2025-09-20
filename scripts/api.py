import os
import joblib
import pandas as pd
from flask import Flask, request, jsonify
from flask_cors import CORS
from google.cloud import storage
import io
import json

app = Flask(__name__)
CORS(app)

# --- Global variables ---
# Test change to verify smart triggers work correctly
model = None
storage_client = None
bucket = None
feature_columns = None

def download_model_and_columns():
    """Downloads the model and feature columns from GCS."""
    global model, feature_columns
    print("Downloading model and feature columns from GCS...")

    # Download model
    model_blob = bucket.blob("model.joblib")
    model_data = model_blob.download_as_string()
    model = joblib.load(io.BytesIO(model_data))
    print("Model downloaded and loaded successfully.")

    # Download feature columns
    columns_blob = bucket.blob("feature_columns.json")
    columns_data = columns_blob.download_as_string()
    feature_columns = json.loads(columns_data)
    print(f"Feature columns loaded. Expecting {len(feature_columns)} features.")


@app.before_request
def load_model_and_columns():
    """Load the model and feature columns before the first request."""
    global storage_client, bucket
    if model is None:
        storage_client = storage.Client()
        bucket = storage_client.get_bucket("zntic-data")
        download_model_and_columns()

@app.route('/predict', methods=['POST'])
def predict():
    """Receives a request and returns a prediction."""
    if not request.json or 'features' not in request.json:
        return jsonify({'error': 'Invalid input: JSON with "features" key is required.'}), 400

    try:
        features = request.json['features']
        
        if len(features) != len(feature_columns):
            return jsonify({'error': f'Invalid number of features. Expected {len(feature_columns)}, got {len(features)}.'}), 400

        # Create a DataFrame with the correct column names
        df = pd.DataFrame([features], columns=feature_columns)

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

@app.route('/info', methods=['GET'])
def info():
    """Returns build and system information."""
    import datetime
    import os

    # Try to get build timestamp from environment or container
    build_date = os.environ.get('BUILD_DATE', 'Unknown')

    # If no build date, try to get from file modification time
    if build_date == 'Unknown':
        try:
            # Get the modification time of this script as fallback
            stat = os.stat(__file__)
            build_date = datetime.datetime.fromtimestamp(stat.st_mtime, tz=datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        except:
            build_date = 'Unknown'

    # Try to get model training date from GCS metadata
    model_training_date = 'Unknown'
    try:
        if bucket:
            model_blob = bucket.blob("model.joblib")
            if model_blob.exists():
                model_blob.reload()
                if model_blob.time_created:
                    model_training_date = model_blob.time_created.strftime('%Y-%m-%d %H:%M:%S UTC')
    except Exception as e:
        print(f"Error getting model training date: {e}")

    return jsonify({
        'api_version': '0.25',
        'build_date': build_date,
        'model_training_date': model_training_date,
        'status': 'operational'
    })

# API endpoints only - web interface served by separate nginx service

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))