# Creates the Kubernetes Service Account (KSA)
resource "kubernetes_service_account" "zntic_gke_sa" {
  metadata {
    name      = "zntic-gke-sa"
    namespace = "default"
  }
}

# Allows the Kubernetes Service Account to impersonate the Google Cloud Service Account
resource "google_service_account_iam_member" "gke_workload_identity_user" {
  service_account_id = "projects/gca-gke-2025/serviceAccounts/zntic-gke-sa@gca-gke-2025.iam.gserviceaccount.com"
  role               = "roles/iam.workloadIdentityUser"
  member             = "serviceAccount:gca-gke-2025.svc.id.goog[default/zntic-gke-sa]"
}
