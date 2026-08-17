FROM nvidia/cuda:12.6.0-base-ubuntu22.04
LABEL org.opencontainers.image.title="nf-core/proteinfold_chai1" \
    org.opencontainers.image.description="Docker image containing all software requirements to run Chai-1 using the nf-core/proteinfold pipeline" \
    org.opencontainers.image.version="1.0.0" \
    org.opencontainers.image.authors="Jose Espinosa-Carrasco <espinosacarrascoj@gmail.com>" \
    org.opencontainers.image.source="https://github.com/nf-core/proteinfold" \
    org.opencontainers.image.licenses="MIT"

# Get latest package listing, install software-properties-common, git and wget.
RUN apt update --quiet \
    && apt install --yes --quiet software-properties-common \
    && apt install --yes --quiet git wget jq \
    && rm -rf /var/lib/apt/lists/*

# Get apt repository of specific Python versions. Then install Python.
# DEBIAN_FRONTEND avoids timezone prompt.
RUN add-apt-repository ppa:deadsnakes/ppa \
    && DEBIAN_FRONTEND=noninteractive apt install --yes --quiet \
        python3.11 \
        python3-pip \
        python3.11-venv \
        python3.11-dev

# Create and activate virtualenv
RUN python3.11 -m venv /chai1_venv
ENV PATH="/chai1_venv/bin:$PATH"

# Upgrade pip inside the venv
RUN pip install --upgrade pip --no-cache-dir

# Clone Chai-1 at a pinned commit for reproducibility
RUN git clone https://github.com/chaidiscovery/chai-lab.git /app/chai1 && \
    cd /app/chai1 && \
    git checkout v0.5.0 && \
    cd -

# Install Chai-1 and its dependencies
RUN pip install --no-cache-dir /app/chai1

# Chai-1 downloads weights at runtime to this path by default;
# expose it as a volume so users can cache weights outside the container.
ENV CHAI_DOWNLOADS_DIR="/data/chai1_weights"
RUN mkdir -p /data/chai1_weights
VOLUME ["/data/chai1_weights"]

# XLA tuning — disable Triton GEMM to avoid long JIT compilation times on some GPUs
ENV XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

# Memory settings suited for A100 80 GB; tune XLA_CLIENT_MEM_FRACTION for smaller GPUs
ENV XLA_PYTHON_CLIENT_PREALLOCATE=true
ENV XLA_CLIENT_MEM_FRACTION=0.95

# Default working directory
WORKDIR /app/chai1