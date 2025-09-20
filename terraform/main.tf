terraform {
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "5.45.2"
    }
    kubernetes = {
      source  = "hashicorp/kubernetes"
      version = "2.38.0"
    }
  }
}

provider "google" {
  project = "gca-gke-2025"
  region  = "us-central1"
}

provider "kubernetes" {
  host                   = "https://${google_container_cluster.primary.endpoint}"
  token                  = data.google_client_config.default.access_token
  cluster_ca_certificate = base64decode(google_container_cluster.primary.master_auth[0].cluster_ca_certificate)
}

data "google_client_config" "default" {}

resource "google_storage_bucket" "data" {
  name     = "zntic-data"
  location = "US"
}

resource "google_storage_bucket" "test_bucket" {
  name     = "zntic-test-bucket"
  location = "US"
}

resource "google_storage_bucket_iam_member" "test_bucket_iam" {
  bucket = google_storage_bucket.test_bucket.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com"
}

resource "google_container_cluster" "primary" {
  name     = "zntic-train"
  location = "us-central1-a"

  # We can't create a cluster with no node pool defined, but we want to only use
  # separately managed node pools. So we create the smallest possible default
  # node pool and immediately delete it.
  remove_default_node_pool = true
  initial_node_count       = 1

  workload_identity_config {
    workload_pool = "gca-gke-2025.svc.id.goog"
  }
}

resource "google_container_node_pool" "default_pool" {
  name       = "default-pool"
  cluster    = google_container_cluster.primary.name
  location   = "us-central1-a"
  node_count = 1

  autoscaling {
    min_node_count = 1
    max_node_count = 3
  }

  node_config {
    machine_type = "e2-medium"
    workload_metadata_config {
      mode = "GKE_METADATA"
    }
  }
}

resource "google_container_node_pool" "gpu_pool" {
  name       = "gpu-pool"
  cluster    = google_container_cluster.primary.name
  location   = "us-central1-a"
  node_count = 0

  autoscaling {
    min_node_count = 0
    max_node_count = 1
  }

  node_config {
    machine_type = "n1-standard-1"
    workload_metadata_config {
      mode = "GKE_METADATA"
    }

    guest_accelerator {
      type  = "nvidia-tesla-t4"
      count = 1
    }
  }
}