FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="Rfantibody" \
    org.opencontainers.image.description="Docker image containing all software requirements to run RFantibody" \
    org.opencontainers.image.version="1.0" \
    org.opencontainers.image.authors="pblas@riojasalud.es" \
    org.opencontainers.image.source="https://github.com/RosettaCommons/RFantibody" \
    org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install --no-install-recommends -y \
    python3.10 python3-pip vim make wget curl \
    && rm -rf /var/lib/apt/lists/*

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"

# Set the working directory to the user's home directory
WORKDIR /home

ENTRYPOINT /bin/bash