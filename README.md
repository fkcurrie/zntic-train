# Zoonotic AI

This project aims to build a cloud-native platform for predicting the zoonotic potential of avian influenza viruses, based on the methods described in the paper "An AI for an AI: identifying zoonotic potential of avian influenza viruses via genomic machine learning."

The platform will be built on Google Cloud Platform (GCP) and is designed to be scalable, reproducible, and automated.

## Features

*   **Prediction API:** A web-based API to submit a viral genome sequence and receive a zoonotic risk score.
*   **Explainable AI:** Identifies and highlights the specific genomic motifs that most influence a prediction.
*   **Automated Training:** A CI/CD pipeline to automatically retrain the model when new data or code is available.
*   **Reproducible Infrastructure:** All cloud infrastructure is defined as code using Terraform, allowing for easy replication.

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