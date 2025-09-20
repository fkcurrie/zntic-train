import os
import joblib
import pandas as pd
from flask import Flask, request, jsonify
from flask_cors import CORS
from google.cloud import storage
from kubernetes import client, config
import io
import json
import logging
import yaml
import uuid

# Configure logging
logging.basicConfig(level=logging.INFO)

app = Flask(__name__)
CORS(app)

# --- Kubernetes Configuration ---
try:
    config.load_incluster_config()
    app.logger.info("Loaded in-cluster Kubernetes config.")
except config.ConfigException:
    try:
        config.load_kube_config()
        app.logger.info("Loaded local kube config.")
    except config.ConfigException:
        app.logger.warning("Could not load any Kubernetes config.")

batch_v1 = client.BatchV1Api()

# --- Model Cache ---
model_cache = {}

def get_model_from_gcs(model_name):
    """Downloads and loads a specific model and its columns from GCS."""
    app.logger.info(f"Cache miss. Attempting to download model '{model_name}' from GCS...")
    try:
        storage_client = storage.Client()
        bucket = storage_client.get_bucket("zntic-data")
        model_path = f"models/{model_name}/model.joblib"
        model_blob = bucket.blob(model_path)
        if not model_blob.exists():
            raise FileNotFoundError(f"Model file not found at {model_path}")
        model_data = model_blob.download_as_string()
        model = joblib.load(io.BytesIO(model_data))
        app.logger.info(f"Model '{model_name}' downloaded and loaded successfully.")
        columns_path = f"models/{model_name}/feature_columns.json"
        columns_blob = bucket.blob(columns_path)
        if not columns_blob.exists():
            raise FileNotFoundError(f"Feature columns not found at {columns_path}")
        columns_data = columns_blob.download_as_string()
        feature_columns = json.loads(columns_data)
        app.logger.info(f"Feature columns for '{model_name}' loaded.")
        return {'model': model, 'feature_columns': feature_columns}
    except Exception as e:
        app.logger.error(f"Failed to download model '{model_name}': {e}", exc_info=True)
        raise

def get_model(model_name):
    """Retrieves a model, using the cache if available."""
    if model_name not in model_cache:
        model_cache[model_name] = get_model_from_gcs(model_name)
    return model_cache[model_name]

@app.route('/', methods=['GET'])
def index():
    """A simple index route to test connectivity."""
    app.logger.info("Root endpoint '/' was hit.")
    return "API is alive!", 200

@app.route('/models', methods=['GET'])
def list_models():
    """Lists available models from the GCS bucket."""
    app.logger.info("Request received for /models endpoint.")
    try:
        storage_client = storage.Client()
        # GCS doesn't have real directories. We list blobs with the prefix
        # and then extract the unique "directory" names.
        blobs = storage_client.list_blobs("zntic-data", prefix='models/')
        
        # Extract the model name from the path: models/MODEL-NAME/file
        model_names = set()
        for blob in blobs:
            # Split the path and check if there are enough parts
            parts = blob.name.split('/')
            if len(parts) > 2:
                model_names.add(parts[1])
        
        sorted_models = sorted(list(model_names))
        app.logger.info(f"Found models: {sorted_models}")
        return jsonify(sorted_models)
    except Exception as e:
        app.logger.error(f"Error listing models: {e}", exc_info=True)
        return jsonify({'error': 'Could not retrieve model list.'}), 500

@app.route('/retrain', methods=['POST'])
def retrain_model():
    """Creates a Kubernetes Job to retrain a model."""
    app.logger.info("Request received for /retrain endpoint.")
    if not request.json:
        return jsonify({'error': 'Invalid input: JSON body is required.'}), 400

    model_type = request.json.get('model_type')
    model_name = request.json.get('model_name')
    hyperparams = request.json.get('hyperparams', {})

    if not model_type or not model_name:
        return jsonify({'error': 'Invalid input: "model_type" and "model_name" are required.'}), 400

    try:
        with open('training-job.yaml', 'r') as f:
            job_manifest = yaml.safe_load(f)

        job_name = f"training-job-{model_name.lower()}-{uuid.uuid4().hex[:6]}"
        job_manifest['metadata']['name'] = job_name
        
        container_args = [
            "--model-type", model_type,
            "--model-name", model_name,
            "--hyperparams", json.dumps(hyperparams)
        ]
        job_manifest['spec']['template']['spec']['containers'][0]['args'] = container_args

        app.logger.info(f"Submitting Kubernetes job '{job_name}' with args: {container_args}")
        batch_v1.create_namespaced_job(body=job_manifest, namespace="default")
        
        return jsonify({
            'status': 'success',
            'message': f"Training job '{job_name}' submitted successfully.",
            'job_name': job_name
        })

    except FileNotFoundError:
        app.logger.error("training-job.yaml not found!")
        return jsonify({'error': 'Internal server error: Training job template not found.'}), 500
    except Exception as e:
        app.logger.error(f"Error creating training job: {e}", exc_info=True)
        return jsonify({'error': f'Failed to create training job: {e}'}), 500

@app.route('/predict', methods=['POST'])
def predict():
    """Receives a request and returns a prediction."""
    if not request.json:
        return jsonify({'error': 'Invalid input: JSON body is required.'}), 400
        
    model_name = request.json.get('model_name')
    features = request.json.get('features')

    if not model_name or not features:
        return jsonify({'error': 'Invalid input: "model_name" and "features" are required.'}), 400

    try:
        model_data = get_model(model_name)
        model = model_data['model']
        feature_columns = model_data['feature_columns']
        
        if len(features) != len(feature_columns):
            return jsonify({'error': f'Invalid number of features. Expected {len(feature_columns)}, got {len(features)}.'}), 400

        df = pd.DataFrame([features], columns=feature_columns)
        prediction = model.predict(df)
        prediction_proba = model.predict_proba(df)

        result = {
            'prediction': 'zoonotic' if prediction[0] == 1 else 'non-zoonotic',
            'confidence': {
                'non-zoonotic': prediction_proba[0][0],
                'zoonotic': prediction_proba[0][1]
            }
        }
        return jsonify(result)

    except FileNotFoundError as e:
        return jsonify({'error': str(e)}), 404
    except Exception as e:
        app.logger.error(f"Prediction error for model '{model_name}': {e}", exc_info=True)
        return jsonify({'error': 'An internal error occurred during prediction.'}), 500

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return "OK", 200

@app.route('/info', methods=['GET'])
def info():
    """Returns build and system information."""
    import datetime
    build_date = os.environ.get('BUILD_DATE', 'Unknown')
    if build_date == 'Unknown':
        try:
            stat = os.stat(__file__)
            build_date = datetime.datetime.fromtimestamp(stat.st_mtime, tz=datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        except Exception as e:
            app.logger.warning(f"Could not determine file modification time: {e}")
            build_date = 'Unknown'
    model_training_date = 'Unknown'
    try:
        info_storage_client = storage.Client()
        info_bucket = info_storage_client.get_bucket("zntic-data")
        model_blob = info_bucket.get_blob("models/test-run-1/model.joblib")
        if model_blob and model_blob.time_created:
            model_training_date = model_blob.time_created.strftime('%Y-%m-%d %H:%M:%S UTC')
        else:
            app.logger.warning("Default model blob not found or has no creation time.")
    except Exception as e:
        app.logger.error(f"Error getting model training date from GCS: {e}", exc_info=True)
        model_training_date = 'Unavailable'
    return jsonify({
        'api_version': '0.41',
        'build_date': build_date,
        'model_training_date': model_training_date,
        'status': 'operational'
    })

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
