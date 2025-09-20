import os
from google.cloud import storage
import uuid

def main():
    """
    A minimal test case to verify GCS access and print the authenticated service account.
    """
    print("✅ Starting GCS test...")

    try:
        # Initialize the GCS client.
        # This will automatically use the service account credentials from the environment.
        storage_client = storage.Client()
        print(f"✅ Authenticated successfully as service account: {storage_client.get_service_account_email()}")
        print(f"✅ Project ID is: {storage_client.project}")

    except Exception as e:
        print(f"❌ Failed to authenticate with GCS: {e}")
        return

    bucket_name = os.environ.get("TEST_BUCKET_NAME")
    if not bucket_name:
        print("❌ TEST_BUCKET_NAME environment variable not set.")
        return

    print(f"\n📝 Starting GCS test on bucket: gs://{bucket_name}")

    try:
        bucket = storage_client.get_bucket(bucket_name)
        file_name = f"test-file-{uuid.uuid4()}.txt"
        blob = bucket.blob(file_name)

        # 1. Create a file
        print(f"  - Creating file: {file_name}...")
        blob.upload_from_string("This is the first line.\n")
        print("  ✅ File created.")

        # 2. Read the file
        print("  - Reading file...")
        content = blob.download_as_string().decode("utf-8")
        print("  ✅ Read content:")
        print(f"---\n{content}---")

        # 3. Append to the file (by re-writing)
        print("  - Appending to file...")
        new_content = content + "This is the second line.\n"
        blob.upload_from_string(new_content)
        print("  ✅ Append complete.")

        # 4. Re-read to verify
        print("  - Re-reading file to verify append...")
        content = blob.download_as_string().decode("utf-8")
        print("  ✅ Read content:")
        print(f"---\n{content}---")

        # 5. Delete the file
        print(f"  - Deleting file...")
        blob.delete()
        print("  ✅ File deleted.")

        # 6. Re-create to ensure deletion worked
        print(f"  - Re-creating file...")
        blob.upload_from_string("This file was re-created.")
        print("  ✅ File re-created.")

        print("\n🎉 GCS test completed successfully!")

    except Exception as e:
        print(f"❌ An error occurred during the GCS test: {e}")

if __name__ == "__main__":
    main()
