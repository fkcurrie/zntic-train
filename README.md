# Zoonotic AI

This project aims to build a cloud-native platform for predicting the zoonotic potential of avian influenza viruses, based on the methods described in the paper "An AI for an AI: identifying zoonotic potential of avian influenza viruses via genomic machine learning."

## Project Configuration & Context

This section provides essential configuration details for any developer or AI assistant joining the project.

*   **GCP Project ID:** `gca-gke-2025`
*   **GKE Cluster:** `zntic-train`
*   **GCS Buckets:**
    *   `zntic-data` (for primary data)
    *   `zntic-test-bucket` (for isolated testing)
*   **CI/CD:** Builds are automatically triggered by pushes to the `main` branch on GitHub via the `zntic-trigger` Cloud Build trigger.
*   **Service Accounts:**
    *   **Google Service Account (GSA):** `zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com`
    *   **Kubernetes Service Account (KSA):** `zntic-gke-sa` (in the `default` namespace)

## Current Status

**✅ GKE Environment Healthy and Configured**

The GKE cluster and associated resources have been rebuilt and are now in a healthy, correctly configured state. The previous GCS permission and autoscaling issues have been resolved.

*   **RESOLVED: GCS Permission Issue:**
    *   **Symptom:** Kubernetes jobs were failing with `403 Forbidden` errors when trying to access GCS.
    *   **Root Cause:** The Google Service Account (GSA) `zntic-gke-sa@...` had the necessary `iam.workloadIdentityUser` role to be impersonated by pods, but the GSA itself **lacked the `roles/storage.admin` permission**. It was authorized to be used, but had no power to access GCS buckets.
    *   **Solution:** The `storage.admin` role was added to the GSA. This fix was embedded directly into the Terraform configuration (`terraform/service-account.tf`) to ensure the permissions are set correctly on infrastructure creation. The `gcs-test-job` now runs successfully.

*   **RESOLVED: Autoscaling Issue:**
    *   **Symptom:** A "Can't scale up nodes" notification was observed in the GCP Console.
    *   **Root Cause:** The `default-pool` node pool was not configured for autoscaling in the Terraform definition (`terraform/main.tf`).
    *   **Solution:** An `autoscaling` block was added to the `default_pool` resource in the Terraform configuration, and the changes were applied to the cluster. The cluster now correctly scales based on workload.

## Cloud Architecture

This project leverages a modern MLOps stack on Google Cloud:

*   **Infrastructure:** [**Terraform**](https://www.terraform.io/) is used to define and manage all cloud resources.
*   **Containerization:** [**Docker**](https://www.docker.com/) is used to package all code.
*   **CI/CD:** [**Google Cloud Build**](https://cloud.google.com/build) automatically builds and pushes container images.
*   **Compute:** [**Google Kubernetes Engine (GKE)**](https://cloud.google.com/kubernetes-engine) orchestrates all containerized workloads, using Workload Identity for secure access to other GCP services.
*   **Data Storage:** [**Google Cloud Storage (GCS)**](https://cloud.google.com/storage) is used for all data and model artifacts.

## Getting Started

1.  **Clone the repository:** `git clone https://github.com/fkcurrie/zntic-train.git`
2.  **Set up authentication:** Ensure you have the `gcloud` CLI installed and authenticated.
3.  **Deploy the infrastructure:** `cd terraform && terraform init && terraform apply`