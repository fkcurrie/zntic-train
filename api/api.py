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
import tensorflow as tf

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

        # Determine model file extension
        model_blob_joblib = bucket.blob(f"models/{model_name}/model.joblib")
        model_blob_keras = bucket.blob(f"models/{model_name}/model.keras")

        if model_blob_keras.exists():
            model_path = model_blob_keras.name
            model_filename = 'model.keras'
            model_blob_keras.download_to_filename(model_filename)
            model = tf.keras.models.load_model(model_filename)
            app.logger.info(f"Keras model '{model_name}' downloaded and loaded successfully.")
        elif model_blob_joblib.exists():
            model_path = model_blob_joblib.name
            model_data = model_blob_joblib.download_as_string()
            model = joblib.load(io.BytesIO(model_data))
            app.logger.info(f"Scikit-learn model '{model_name}' downloaded and loaded successfully.")
        else:
            raise FileNotFoundError(f"No model file (.joblib or .keras) found for model '{model_name}'")

        columns_path = f"models/{model_name}/feature_columns.json"
        columns_blob = bucket.blob(columns_path)
        if not columns_blob.exists():
            raise FileNotFoundError(f"Feature columns not found at {columns_path}")
        columns_data = columns_blob.download_as_string()
        feature_columns = json.loads(columns_data)
        app.logger.info(f"Feature columns for '{model_name}' loaded.")
        
        return {'model': model, 'feature_columns': feature_columns, 'type': 'keras' if model_blob_keras.exists() else 'sklearn'}
    except Exception as e:
        app.logger.error(f"Failed to download model '{model_name}': {e}", exc_info=True)
        raise

def get_model(model_name):
    """Retrieves a model, using the cache if available."""
    if model_name not in model_cache:
        model_cache[model_name] = get_model_from_gcs(model_name)
    return model_cache[model_name]

@app.route('/models', methods=['GET'])
def list_models():
    """Lists available models from the GCS bucket."""
    app.logger.info("Request received for /models endpoint.")
    try:
        storage_client = storage.Client()
        blobs = storage_client.list_blobs("zntic-data", prefix='models/')
        model_names = set()
        for blob in blobs:
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
        # Choose the correct job template based on model type
        if model_type == 'NeuralNetwork':
            template_file = 'training-job-gpu.yaml'
        else:
            template_file = 'training-job.yaml'
            
        with open(template_file, 'r') as f:
            job_manifest = yaml.safe_load(f)

        job_name = f"training-job-{model_name.lower().replace('_', '-')}-{uuid.uuid4().hex[:6]}"
        job_manifest['metadata']['name'] = job_name
        
        container_args = [
            "--model-type", model_type,
            "--model-name", model_name,
            "--hyperparams", json.dumps(hyperparams)
        ]
        job_manifest['spec']['template']['spec']['containers'][0]['args'] = container_args

        app.logger.info(f"Submitting Kubernetes job '{job_name}' with args: {container_args} using template {template_file}")
        batch_v1.create_namespaced_job(body=job_manifest, namespace="default")
        
        return jsonify({
            'status': 'success',
            'message': f"Training job '{job_name}' submitted successfully.",
            'job_name': job_name
        })

    except FileNotFoundError:
        app.logger.error(f"{template_file} not found!")
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
        model_type = model_data['type']
        
        if len(features) != len(feature_columns):
            return jsonify({'error': f'Invalid number of features. Expected {len(feature_columns)}, got {len(features)}.'}), 400

        df = pd.DataFrame([features], columns=feature_columns)
        
        if model_type == 'keras':
            prediction_proba = model.predict(df)
            prediction = (prediction_proba > 0.5).astype(int)
            # Keras returns probabilities in a different shape
            prediction_proba_formatted = [[1 - prediction_proba[0][0], prediction_proba[0][0]]]
        else: # sklearn
            prediction = model.predict(df)
            prediction_proba_formatted = model.predict_proba(df)

        result = {
            'prediction': 'zoonotic' if prediction[0] == 1 else 'non-zoonotic',
            'confidence': {
                'non-zoonotic': prediction_proba_formatted[0][0],
                'zoonotic': prediction_proba_formatted[0][1]
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

# Other endpoints like /info are omitted for brevity but would be here

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))# Trivial change to trigger build
