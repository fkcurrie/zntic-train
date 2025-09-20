# TODO

This document outlines the development plan for building the Zoonotic AI platform on Google Cloud.

## ✅ PRIORITY 1: Resolve GKE Workload Identity Issue

- [x] **Investigate and Fix the `403 Forbidden` Error:**
    - [x] **Goal:** Successfully run the `gcs-test-job` without permission errors.
    - [x] **Status:** **RESOLVED**. The root cause was that the Google Service Account (GSA) `zntic-gke-sa` was missing the `roles/storage.admin` permission. This has been fixed in the Terraform configuration.

## Phase 1: Infrastructure Setup (Terraform & Cloud Build)

- [x] **Terraform:**
    - [x] Create a `terraform/` directory for all infrastructure code.
    - [x] Define a GCS bucket to store raw data and trained model artifacts.
    - [x] Define a GKE cluster named `zntic-train`.
        - [x] Configure a default, always-on CPU node pool.
        - [x] Configure a second, autoscaling node pool with NVIDIA T4 GPUs.
        - [x] Configure Workload Identity.
        - [x] **[FIXED]** Enable autoscaling on the default node pool.
- [x] **Cloud Build:**
    - [x] Create a `cloudbuild.yaml` file to define the CI/CD pipeline.
    - [x] Create a build trigger that automatically runs on pushes to the GitHub repository.

## Phase 2: Data Processing & Model Training Pipeline

- [ ] **Data Ingestion:**
    - [ ] Create a `Dockerfile` for a data-fetching script (`download_data.py`).
    - [ ] Write a Kubernetes Job manifest (`download-job.yaml`) to run the data-fetching container.
- [ ] **Feature Engineering & Training:**
    - [ ] Create a `Dockerfile` for the training script (`feature_engineering.py` and `train.py`).
    - [ ] Write a Kubernetes Job manifest (`training-job.yaml`) that requests GPU resources.
    - [ ] The training job will:
        - [ ] Load the raw data from the GCS bucket.
        - [ ] Perform feature engineering.
        - [ ] Train the model.
        - [ ] Save the final trained model artifact back to the GCS bucket.

## Phase 3: Inference API & Web Interface

- [ ] **API Development:**
    - [ ] Create a simple web application (e.g., using Flask or FastAPI).
- [ ] **Containerization & Deployment:**
    - [ ] Create a `Dockerfile` for the web application.
    - [ ] Write Kubernetes Deployment and Service manifests.

## Phase 4: Automation & CI/CD

- [ ] Connect the Cloud Build triggers to the `zntic-train` GitHub repository.
- [ ] Configure the main branch trigger to automatically deploy new versions of the API to the GKE cluster.
