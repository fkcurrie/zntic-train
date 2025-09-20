# Zoonotic AI

This project aims to build a cloud-native platform for predicting the zoonotic potential of avian influenza viruses, based on the methods described in the paper "An AI for an AI: identifying zoonotic potential of avian influenza viruses via genomic machine learning."

The platform is built on Google Cloud Platform (GCP) and is designed to be scalable, reproducible, and automated.

## Project Status

**Infrastructure:** The core infrastructure, consisting of a GKE cluster and a GCS bucket, has been successfully provisioned using Terraform.

**CI/CD:** A Cloud Build trigger is in place to automatically build and push container images to GCR when changes are pushed to the `main` branch.

**Data Ingestion:** A containerized data download script has been developed. The process of running this as a Kubernetes job on the GKE cluster to populate the GCS bucket is currently in progress.

## Cloud Architecture

This project leverages a modern MLOps stack on Google Cloud:

*   **Infrastructure:** [**Terraform**](https://www.terraform.io/) is used to define and manage all cloud resources, including the GKE cluster and GCS buckets. This ensures the entire environment can be created and destroyed repeatably.
*   **Containerization:** [**Docker**](https://www.docker.com/) is used to package the data processing, model training, and web application code into portable containers.
*   **CI/CD:** [**Google Cloud Build**](https://cloud.google.com/build) automatically builds and tests the container images whenever changes are pushed to the GitHub repository.
*   **Compute for Training & Serving:** [**Google Kubernetes Engine (GKE)**](https://cloud.google.com/kubernetes-engine) orchestrates all containerized workloads. The cluster is configured with:
    *   A standard CPU-based node pool for running the 24/7 web application.
    *   An autoscaling GPU (NVIDIA T4) node pool that scales down to zero to run computationally intensive training jobs cost-effectively.
*   **Data Storage:** [**Google Cloud Storage (GCS)**](https://cloud.google.com/storage) provides a scalable and durable location for storing raw genomic data, processed datasets, and trained model artifacts.

## Technology Stack

*   **Machine Learning:** Python, TensorFlow/Keras, Scikit-learn, Pandas, Biopython
*   **Cloud & DevOps:** Google Cloud Platform (GKE, Cloud Build, GCS), Terraform, Docker, GitHub

## Getting Started

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/fkcurrie/zntic-train.git
    cd zntic-train
    ```

2.  **Set up authentication:**
    Ensure you have the `gcloud` CLI installed and authenticated with your GCP account.

3.  **Deploy the infrastructure:**
    Navigate to the `terraform` directory and run the following commands:
    ```bash
    terraform init
    terraform apply
    ```
