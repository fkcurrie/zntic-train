# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a cloud-native MLOps platform for predicting zoonotic potential of avian influenza viruses using genomic machine learning. The project runs entirely on Google Cloud Platform with containerized workloads orchestrated by Kubernetes.

## Essential Configuration

- **GCP Project:** `gca-gke-2025`
- **GKE Cluster:** `zntic-train`
- **Primary GCS Bucket:** `zntic-data`
- **Test GCS Bucket:** `zntic-test-bucket`
- **Google Service Account:** `zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com`
- **Kubernetes Service Account:** `zntic-gke-sa` (default namespace)
- **Container Registry:** `gcr.io/gca-gke-2025/`

## Architecture

The platform follows a microservices architecture with these main components:

1. **Data Pipeline** (`scripts/download_data.py`) - Fetches and prepares genomic data
2. **Feature Engineering** (`scripts/feature_engineering.py`) - Processes genomic sequences  
3. **Model Training** (`scripts/train.py`) - Trains ML models on processed features
4. **Inference API** (`scripts/api.py`) - Flask REST API for model predictions
5. **Infrastructure** (`terraform/`) - Terraform configs for GCP resources

All components run as containerized Kubernetes Jobs or Deployments, with data and models stored in GCS.

## Common Commands

### Infrastructure Management
```bash
# Deploy/update infrastructure
cd terraform && terraform init && terraform apply

# Destroy infrastructure (careful!)
cd terraform && terraform destroy
```

### Container Operations
```bash
# Build all containers (via Cloud Build)
gcloud builds submit --config=cloudbuild.yaml

# Test GCS connectivity
kubectl apply -f gcs-test-job.yaml
kubectl logs -f job/gcs-test-job
```

### ML Pipeline
```bash
# Run data download
kubectl apply -f download-job.yaml

# Run training pipeline  
kubectl apply -f training-job.yaml

# Deploy inference API
kubectl apply -f api-deployment.yaml
kubectl apply -f api-service.yaml
```

### Development Workflow
```bash
# Install Python dependencies
pip install -r requirements.txt

# Local testing of scripts
python scripts/gcs_test.py  # Test GCS access
python scripts/api.py       # Run API locally
```

## Key Technical Details

- **Workload Identity** is configured for secure GCS access from Kubernetes pods
- **Autoscaling** is enabled on the default node pool
- **GPU node pool** with T4 GPUs available for training workloads
- **Cloud Build** automatically triggers on main branch pushes
- **Feature columns** are saved/loaded via `feature_columns.json` in GCS for API consistency
- **Flask API** expects POST requests to `/predict` with genomic features

## Container Images

The build process creates four main images:
- `gcr.io/gca-gke-2025/download-data` - Data fetching
- `gcr.io/gca-gke-2025/training` - Model training  
- `gcr.io/gca-gke-2025/zoonotic-api` - Inference API
- `gcr.io/gca-gke-2025/gcs-test` - GCS connectivity testing

## Dependencies

Core Python packages (see `requirements.txt`):
- `biopython` - Genomic sequence processing
- `pandas`, `scikit-learn` - Data processing and ML
- `google-cloud-storage` - GCS integration
- `Flask`, `gunicorn` - Web API framework