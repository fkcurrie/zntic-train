# TODO

This document outlines the development plan for building the Zoonotic AI platform on Google Cloud.

## Phase 1: Infrastructure Setup (Terraform & Cloud Build)

- [ ] **Terraform:**
    - [ ] Create a `terraform/` directory for all infrastructure code.
    - [ ] Define a GCS bucket to store raw data and trained model artifacts.
    - [ ] Define a GKE cluster named `zntic-train`.
        - [ ] Configure a default, always-on CPU node pool for the web application.
        - [ ] Configure a second, autoscaling node pool with NVIDIA T4 GPUs for training, with a minimum size of 0.
- [ ] **Cloud Build:**
    - [ ] Create a `cloudbuild.yaml` file to define the CI/CD pipeline.
    - [ ] Create a build trigger that automatically runs on pushes to the GitHub repository.

## Phase 2: Data Processing & Model Training Pipeline

- [ ] **Data Ingestion:**
    - [ ] Create a `Dockerfile` for a data-fetching script (`download_data.py`).
    - [ ] Write a Kubernetes Job manifest (`download-job.yaml`) to run the data-fetching container on the GKE cluster. This job will download the initial dataset and save it to the GCS bucket.
- [ ] **Feature Engineering & Training:**
    - [ ] Create a `Dockerfile` for the training script (`feature_engineering.py` and `train.py`). This container will have all ML dependencies.
    - [ ] Write a Kubernetes Job manifest (`training-job.yaml`) that requests GPU resources from the cluster.
    - [ ] The training job will:
        - [ ] Load the raw data from the GCS bucket.
        - [ ] Perform feature engineering.
        - [ ] Train the model using the GPU node pool.
        - [ ] Save the final trained model artifact back to the GCS bucket.

## Phase 3: Inference API & Web Interface

- [ ] **API Development:**
    - [ ] Create a simple web application (e.g., using Flask or FastAPI) that:
        - [ ] Loads the trained model from the GCS bucket at startup.
        - [ ] Provides a `/predict` endpoint that accepts a FASTA sequence.
        - [ ] Returns a JSON response with the zoonotic risk score and influential motifs.
- [ ] **Containerization & Deployment:**
    - [ ] Create a `Dockerfile` for the web application.
    - [ ] Update `cloudbuild.yaml` to build and push this container image to Google Artifact Registry.
    - [ ] Write a Kubernetes Deployment manifest (`api-deployment.yaml`) to run the web application on the default CPU node pool.
    - [ ] Write a Kubernetes Service manifest (`api-service.yaml`) to expose the web application to the internet via a LoadBalancer.

## Phase 4: Automation & CI/CD

- [ ] Connect the Cloud Build triggers to the `zntic-train` GitHub repository.
- [ ] Configure the main branch trigger to automatically deploy new versions of the API to the GKE cluster.
