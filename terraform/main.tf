terraform {
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 5.0"
    }
  }
}

provider "google" {
  project = "gca-gke-2025"
  region  = "us-central1"
}

# A GCS bucket to store the training data
resource "google_storage_bucket" "training_data" {
  name          = "zntic-train-data"
  location      = "US-CENTRAL1"
  force_destroy = true # Set to true for easier cleanup in a dev environment

  uniform_bucket_level_access = true
}

resource "google_storage_bucket_iam_member" "gke_storage_access" {
  bucket = google_storage_bucket.training_data.name
  role   = "roles/storage.admin"
  member = "serviceAccount:zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com"
}

resource "google_project_iam_member" "gke_artifact_reader" {
  project = "gca-gke-2025"
  role    = "roles/artifactregistry.reader"
  member  = "serviceAccount:zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com"
}

resource "google_container_cluster" "zntic_train" {
  name     = "zntic-train"
  location = "us-central1"

  # We can't remove the default node pool, but we can make it small.
  remove_default_node_pool = true
  initial_node_count       = 1

  # Enable Workload Identity
  workload_identity_config {
    workload_pool = "gca-gke-2025.svc.id.goog"
  }
}

# Default node pool for general workloads (non-GPU)
resource "google_container_node_pool" "default_pool" {
  name       = "default-cpu-pool"
  cluster    = google_container_cluster.zntic_train.name
  location   = google_container_cluster.zntic_train.location
  node_count = 1

  node_config {
    machine_type = "e2-medium"
    service_account = "zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com"
    oauth_scopes = [
      "https://www.googleapis.com/auth/cloud-platform"
    ]
  }
}

# Autoscaling node pool with T4 GPUs for training jobs
resource "google_container_node_pool" "gpu_training_pool" {
  name     = "t4-gpu-training-pool"
  cluster  = google_container_cluster.zntic_train.name
  location = google_container_cluster.zntic_train.location

  autoscaling {
    min_node_count = 0
    max_node_count = 3 # Can be adjusted based on workload
  }
  
  management {
    auto_repair  = true
    auto_upgrade = true
  }

  node_config {
    machine_type = "n1-standard-4" # A good general-purpose machine type for a T4
    service_account = "zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com"
    
    guest_accelerator {
      type  = "nvidia-tesla-t4"
      count = 1
    }

    # Taint nodes to only schedule pods that explicitly request GPUs
    taint {
      key    = "nvidia.com/gpu"
      value  = "present"
      effect = "NO_SCHEDULE"
    }

    oauth_scopes = [
      "https://www.googleapis.com/auth/cloud-platform"
    ]
  }
}
