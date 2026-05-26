FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --timeout 100 --retries 10 pandas biopython

CMD ["python"]
