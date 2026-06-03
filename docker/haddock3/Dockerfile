FROM python:3.11

RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir haddock3

WORKDIR /data

ENTRYPOINT ["haddock3"]