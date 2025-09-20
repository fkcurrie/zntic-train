import os
from google.cloud import storage
import google.auth

def gcs_test():
    """
    Performs a series of GCS operations to test permissions.
    """
    try:
        # --- 1. AUTHENTICATION CHECK ---
        credentials, project_id = google.auth.default()
        if hasattr(credentials, 'service_account_email'):
            print(f"✅ Authenticated successfully as service account: {credentials.service_account_email}")
        else:
            print(f"✅ Authenticated successfully with user credentials.")
        print(f"✅ Project ID is: {project_id}")

        # --- 2. SETUP ---
        bucket_name = os.environ.get("TEST_BUCKET_NAME")
        if not bucket_name:
            print("❌ ERROR: TEST_BUCKET_NAME environment variable not set.")
            return

        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        file_name = "test-file.txt"
        blob = bucket.blob(file_name)

        print(f"\n📝 Starting GCS test on bucket: gs://{bucket_name}")

        # --- 3. CREATE ---
        print(f"  - Creating file: {file_name}...")
        blob.upload_from_string("This is the first line.\n")
        print(f"  ✅ File created.")

        # --- 4. READ ---
        print(f"  - Reading file...")
        content = blob.download_as_string().decode("utf-8")
        print(f"  ✅ Read content:\n---\n{content.strip()}\n---")

        # --- 5. APPEND ---
        print(f"  - Appending to file...")
        new_content = content + "This is the second line.\n"
        blob.upload_from_string(new_content)
        print(f"  ✅ Append complete.")

        # --- 6. RE-READ ---
        print(f"  - Re-reading file to verify append...")
        content = blob.download_as_string().decode("utf-8")
        print(f"  ✅ Read content:\n---\n{content.strip()}\n---")
        assert "second line" in content

        # --- 7. DELETE ---
        print(f"  - Deleting file...")
        blob.delete()
        print(f"  ✅ File deleted.")

        # --- 8. RE-CREATE ---
        print(f"  - Re-creating file...")
        blob.upload_from_string("This is the final line.\n")
        print(f"  ✅ File re-created.")

        print("\n🎉 GCS test completed successfully!")

    except Exception as e:
        print(f"\n❌ An error occurred during the GCS test: {e}")
        # Re-raise the exception to ensure the Kubernetes job fails.
        raise

if __name__ == "__main__":
    gcs_test()
