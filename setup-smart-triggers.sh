#!/bin/bash

# Smart Cloud Build Triggers Setup
# Only rebuild services when their specific files change

echo "Setting up smart Cloud Build triggers with path filters..."

# API trigger - only rebuild when API-related files change
echo "Creating API auto-trigger..."
gcloud builds triggers create github \
  --repo-owner=fkcurrie \
  --repo-name=zntic-train \
  --branch-pattern="^main$" \
  --build-config=cloudbuild-api.yaml \
  --included-files="scripts/api.py,api.Dockerfile,requirements.txt" \
  --name=api-auto-trigger \
  --description="Auto-build API service when API files change" \
  --region=us-central1

# Training trigger - only rebuild when training-related files change
echo "Creating Training auto-trigger..."
gcloud builds triggers create github \
  --repo-owner=fkcurrie \
  --repo-name=zntic-train \
  --branch-pattern="^main$" \
  --build-config=cloudbuild-training.yaml \
  --included-files="scripts/train.py,scripts/feature_engineering.py,training.Dockerfile,requirements.txt" \
  --name=training-auto-trigger \
  --description="Auto-build training service when training files change" \
  --region=us-central1

# Download trigger - only rebuild when download-related files change
echo "Creating Download auto-trigger..."
gcloud builds triggers create github \
  --repo-owner=fkcurrie \
  --repo-name=zntic-train \
  --branch-pattern="^main$" \
  --build-config=cloudbuild-download.yaml \
  --included-files="scripts/download_data.py,scripts/Dockerfile.download,requirements.txt" \
  --name=download-auto-trigger \
  --description="Auto-build download service when download files change" \
  --region=us-central1

echo "Smart triggers setup complete!"
echo ""
echo "Now when you push changes:"
echo "- Only API will rebuild if you change scripts/api.py or api.Dockerfile"
echo "- Only training will rebuild if you change scripts/train.py or training files"
echo "- Only download will rebuild if you change scripts/download_data.py"
echo "- If you change requirements.txt, all services will rebuild (as expected)"