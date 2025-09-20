FROM python:3.9-slim

WORKDIR /app

COPY scripts/download_data.py .
COPY requirements.txt .

RUN apt-get update && apt-get install -y curl && \
    curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' && \
    curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat' && \
    chmod +x datasets dataformat && \
    mv datasets dataformat /usr/local/bin/ && \
    pip install --no-cache-dir -r requirements.txt

CMD ["python", "download_data.py"]