#!/bin/bash

# Complete Cloud Build Triggers Setup
# Run this script AFTER completing GitHub OAuth authorization

echo "=== Complete Cloud Build Triggers Setup ==="
echo "This script will create the repository connection and smart triggers."
echo ""

# Check if GitHub connection is properly authorized
echo "Checking GitHub connection status..."
CONNECTION_STATUS=$(gcloud builds connections describe zntic-github-connection --region=us-central1 --format="value(installationState.stage)" 2>/dev/null)

if [ "$CONNECTION_STATUS" != "COMPLETE" ]; then
    echo "❌ ERROR: GitHub connection is not properly authorized."
    echo "Current status: $CONNECTION_STATUS"
    echo ""
    echo "Please complete GitHub OAuth authorization first:"
    echo "Visit: https://accounts.google.com/AccountChooser?continue=https%3A%2F%2Fconsole.cloud.google.com%2Fm%2Fgcb%2Fgithub%2Flocations%2Fus-central1%2Foauth_v2%3Fconnection_name%3Dprojects%252F764460891170%252Flocations%252Fus-central1%252Fconnections%252Fzntic-github-connection"
    echo ""
    echo "After authorization, run this script again."
    exit 1
fi

echo "✅ GitHub connection is authorized!"
echo ""

# Create repository connection
echo "Creating repository connection..."
gcloud builds repositories create zntic-train \
  --connection=zntic-github-connection \
  --region=us-central1 \
  --remote-uri=https://github.com/fkcurrie/zntic-train.git

if [ $? -ne 0 ]; then
    echo "❌ Failed to create repository connection"
    exit 1
fi

echo "✅ Repository connection created!"
echo ""

# Create smart triggers with path filters
echo "Creating smart Cloud Build triggers with path filters..."
echo ""

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

# Removed GCS test trigger as no longer needed

echo ""
echo "🎉 Smart triggers setup complete!"
echo ""
echo "Now when you push changes to the main branch:"
echo "✅ Only API will rebuild if you change scripts/api.py or api.Dockerfile"
echo "✅ Only training will rebuild if you change scripts/train.py or training files"
echo "✅ Only download will rebuild if you change scripts/download_data.py"
echo "✅ If you change requirements.txt, all services will rebuild (as expected)"
echo ""
echo "You can verify the triggers in the Cloud Console:"
echo "https://console.cloud.google.com/cloud-build/triggers?project=gca-gke-2025"