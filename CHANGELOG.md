# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - 2025-09-20

### Fixed

- Resolved all issues with the data ingestion pipeline. The `download-data` job now runs successfully on GKE.

### Changed

- Updated `download-job.yaml` to correctly mount service account credentials.
- Updated `download_data.py` to explicitly use service account credentials.
- Updated `requirements.txt` to include all necessary dependencies.
- Updated `Dockerfile` to correctly install `ncbi-datasets-cli`.

## [0.3.0] - 2025-09-20

### Added

- `training.Dockerfile` for containerizing the training pipeline.
- `feature_engineering.py` and `train.py` scripts.
- `training-job.yaml` to run the training pipeline on GKE.

### Changed

- Updated `cloudbuild.yaml` to build and push the training image.
- Updated `requirements.txt` with new dependencies.

## [0.2.0] - 2025-09-19

### Added

- Terraform configuration for GCS bucket and GKE cluster.
- `cloudbuild.yaml` for CI/CD pipeline.
- `Dockerfile` for data ingestion.
- Kubernetes Job manifest for downloading data.

## [0.1.0] - 2025-09-19

### Added

- Initial project setup.
- README.md, TODO.md, and CHANGELOG.md files.