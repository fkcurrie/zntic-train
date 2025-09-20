# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file into the container at /app
COPY requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the scripts into the container
COPY scripts/feature_engineering.py /app/
COPY scripts/train.py /app/

# Make the scripts executable
RUN chmod +x /app/feature_engineering.py
RUN chmod +x /app/train.py