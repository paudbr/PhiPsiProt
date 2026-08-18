FROM python:3.11-slim

RUN pip install --no-cache-dir \
    --timeout 120 \
    -i https://pypi.org/simple \
    biopython

WORKDIR /app
