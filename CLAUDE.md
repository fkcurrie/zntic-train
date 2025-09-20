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

1. **Data Pipeline** (`scripts/download_data.py`) - Fetches and prepares genomic data from GISAID
2. **Feature Engineering** (`scripts/feature_engineering.py`) - Extracts 64D nucleotide triplet features from genomic sequences
3. **Model Training** (`scripts/train.py`) - Trains Random Forest classifier on processed features
4. **Inference API** (`scripts/api.py`) - Flask REST API for model predictions (`/predict`, `/health`)
5. **Web Interface** (`web/`) - Academic research interface for genomic analysis
6. **Infrastructure** (`terraform/`) - Terraform configs for GCP resources

### Current Deployment Architecture

- **API Service**: `http://34.30.222.113` - ML inference endpoints
- **Web Interface**: `http://34.42.26.58` - Academic research portal
- **Automatic Deployment**: Cloud Build triggers for code changes
- **Container Orchestration**: Kubernetes with rolling deployments
- **Data Storage**: GCS buckets for models, datasets, and feature metadata

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
# Build all containers individually (via Cloud Build)
gcloud builds submit --config=cloudbuild-api.yaml      # ML inference API
gcloud builds submit --config=cloudbuild-web.yaml      # Web interface
gcloud builds submit --config=cloudbuild-training.yaml # Training pipeline
gcloud builds submit --config=cloudbuild-download.yaml # Data pipeline

# Build all containers (legacy monolithic config)
gcloud builds submit --config=cloudbuild.yaml

# Test GCS connectivity
kubectl apply -f gcs-test-job.yaml
kubectl logs -f job/gcs-test-job

# Check auto-trigger status
gcloud builds triggers list
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

# Deploy web interface
kubectl apply -f web-deployment.yaml

# Check service status
kubectl get services
kubectl get deployments
kubectl rollout status deployment/zoonotic-api-deployment
kubectl rollout status deployment/zoonotic-web-deployment
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
- **Cloud Build** automatically triggers on main branch pushes via configured triggers
- **Feature columns** are saved/loaded via `feature_columns.json` in GCS for API consistency
- **Flask API** expects POST requests to `/predict` with 64-dimensional genomic features
- **Web Interface** provides academic research portal with sample data and automatic analysis
- **Microservices** architecture separates ML API from web interface for scalability
- **Auto-deployment** triggers rebuild web interface when files in `web/` directory change

## Container Images

The build process creates five main images:
- `gcr.io/gca-gke-2025/download-data` - Data fetching from GISAID
- `gcr.io/gca-gke-2025/training` - Model training (Random Forest)
- `gcr.io/gca-gke-2025/zoonotic-api` - Flask ML inference API
- `gcr.io/gca-gke-2025/zoonotic-web` - Nginx-based web interface
- `gcr.io/gca-gke-2025/gcs-test` - GCS connectivity testing

## Auto-Deployment Triggers

Cloud Build triggers automatically deploy services when code changes:
- **Web Interface**: Triggers on changes to `web/` files (index.html, script.js, styles.css, nginx.conf, web.Dockerfile)
- **API Service**: Triggers on changes to API-related files
- **Training Pipeline**: Triggers on changes to training scripts

Trigger management:
```bash
# List all triggers
gcloud builds triggers list

# Import trigger from YAML
gcloud builds triggers import --source=web-trigger.yaml
```

## Dependencies

Core Python packages (see `requirements.txt`):
- `biopython` - Genomic sequence processing
- `pandas`, `scikit-learn` - Data processing and ML
- `google-cloud-storage` - GCS integration
- `Flask`, `gunicorn` - Web API framework
- `flask-cors` - Cross-origin resource sharing for web interface

## Current Status (v0.25)

### Completed Features
- ✅ Microservices architecture with separate API and web services
- ✅ Academic-styled web interface for research use
- ✅ Automatic deployment triggers for continuous integration
- ✅ Random Forest model for zoonotic risk prediction
- ✅ 64-dimensional nucleotide triplet feature engineering
- ✅ GCS-based model and data storage
- ✅ Kubernetes orchestration with rolling deployments
- ✅ Sample data generation for testing and validation

### Known Issues
- ⚠️ Single model architecture (Random Forest only)
- ⚠️ Fixed feature engineering approach (triplet nucleotides)
- ⚠️ No model retraining capability through web interface
- ⚠️ Limited to batch predictions (no real-time streaming)

## Future Enhancements (Roadmap)

### Model Enhancement (Priority 1)
- 🎯 **Multi-Model Architecture**: Support for multiple ML algorithms (XGBoost, Neural Networks, SVM)
- 🎯 **Dynamic Feature Engineering**: Configurable feature extraction methods
- 🎯 **Model Comparison**: A/B testing framework for model performance
- 🎯 **Online Learning**: Incremental model updates with new data
- 🎯 **Model Versioning**: Track and rollback model deployments

### Advanced Feature Engineering (Priority 1)
- 🔬 **Sequence-based Features**: K-mer analysis, GC content, codon usage bias
- 🔬 **Protein-level Features**: Amino acid composition, secondary structure prediction
- 🔬 **Phylogenetic Features**: Evolutionary distance metrics
- 🔬 **Metadata Integration**: Geographic, temporal, and host species features
- 🔬 **Deep Learning Features**: Learned representations from sequence autoencoders

### Web Interface Enhancements (Priority 2)
- 🖥️ **Model Retraining Interface**: Upload datasets and retrain models via web UI
- 🖥️ **Feature Engineering Pipeline**: Configure feature extraction through web interface
- 🖥️ **Batch Processing**: Upload and process multiple sequences
- 🖥️ **Results Visualization**: Interactive charts, ROC curves, feature importance plots
- 🖥️ **Export Capabilities**: Download results in various formats (CSV, PDF, JSON)

### Data Pipeline Improvements (Priority 2)
- 📊 **Real-time Data Ingestion**: Stream processing for new sequences
- 📊 **Data Quality Monitoring**: Automated data validation and quality checks
- 📊 **Federated Learning**: Collaborate with multiple data sources
- 📊 **Data Augmentation**: Synthetic sequence generation for training
- 📊 **Version Control**: Track dataset versions and lineage

### Infrastructure & Operations (Priority 3)
- 🚀 **Multi-region Deployment**: Global availability and disaster recovery
- 🚀 **API Rate Limiting**: Throttling and quota management
- 🚀 **Model Serving Optimization**: TensorFlow Serving, ONNX Runtime
- 🚀 **Monitoring & Alerting**: Comprehensive observability stack
- 🚀 **Security Hardening**: Authentication, authorization, audit logging