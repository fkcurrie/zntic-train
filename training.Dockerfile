FROM python:3.9-slim

WORKDIR /app

COPY scripts/feature_engineering.py .
COPY scripts/train.py .
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "train.py"]
