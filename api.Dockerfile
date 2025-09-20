# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file into the container at /app
COPY requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application scripts into the container
COPY scripts/api.py /app/
COPY training-job.yaml /app/

# Copy the web files (though not served by this container, they might be needed for context)
COPY web/ /app/web/

# Expose the port the app runs on
EXPOSE 8080

# Run the application using gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "api:app"]