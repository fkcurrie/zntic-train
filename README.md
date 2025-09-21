# GISAID Assessment Tool

Welcome to the GISAID Assessment Tool, a cloud-native platform for predicting the zoonotic potential of avian influenza viruses. This tool provides a user-friendly interface to run predictions, manage a library of machine learning models, and even train new models on demand.

## Features

This platform is a complete MLOps solution that allows you to:

*   **Predict Zoonotic Risk:** Input genomic feature vectors and get a probabilistic risk assessment from a variety of machine learning models.
*   **Select Your Model:** Choose from a library of pre-trained models, including Logistic Regression, Random Forest, SVC, Gradient Boosting, and Neural Networks.
*   **Train New Models:** Use the Retraining Dashboard to configure and launch new training jobs on the Kubernetes cluster. Experiment with different model types and hyperparameters.
*   **GPU-Accelerated Training:** Neural Network models are automatically trained on a dedicated GPU-enabled node pool for maximum efficiency.
*   **Automated CI/CD:** Changes to the source code automatically trigger builds and deployments for the relevant components.

## How to Use the Tool

The primary way to interact with this platform is through the web interface.

1.  **Prediction Interface:**
    *   Select a pre-trained model from the dropdown menu.
    *   Input your 64-dimensional genomic feature vector.
    *   Click "Run Analysis" to see the prediction and confidence score.

2.  **Retraining Dashboard:**
    *   Navigate to the "Retrain Models" page.
    *   Give your new model a unique name.
    *   Select a model type (e.g., `RandomForestClassifier`, `NeuralNetwork`).
    *   Adjust the default hyperparameters in the JSON editor as needed.
    *   Click "Start Training Job" to launch the training process on the GKE cluster. The new model will appear in the prediction dropdown once training is complete.

## Cloud Architecture

This project leverages a modern MLOps stack on Google Cloud for scalability and automation:

*   **Infrastructure as Code:** [**Terraform**](https://www.terraform.io/) defines and manages all cloud resources, including the GKE cluster and GCS buckets.
*   **Containerization:** [**Docker**](https://www.docker.com/) packages the API, training, and web components into reproducible container images.
*   **CI/CD:** [**Google Cloud Build**](https://cloud.google.com/build) uses a "smart trigger" system. Pushes to the repository only build the specific components that have changed.
*   **Compute:** [**Google Kubernetes Engine (GKE)**](https://cloud.google.com/kubernetes-engine) orchestrates all containerized workloads, with separate node pools for standard and GPU-intensive (NVIDIA T4) tasks.
*   **Secure Access:** GKE Workload Identity provides pods with secure, keyless access to other GCP services.
*   **Data & Model Storage:** [**Google Cloud Storage (GCS)**](https://cloud.google.com/storage) serves as the central repository for all data, model artifacts, and build logs.

---
*This README was last updated by your AI assistant.*
# Trivial change
# Final test
