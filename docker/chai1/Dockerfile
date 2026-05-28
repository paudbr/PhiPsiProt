FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

LABEL org.opencontainers.image.title="proteinfold_chai1" \
    org.opencontainers.image.description="Docker image containing all software requirements to run Chai-1" \
    org.opencontainers.image.version="1.3.0" \
    org.opencontainers.image.authors="pblas@riojasalud.es" \
    org.opencontainers.image.source="https://github.com/nf-core/proteinfold" \
    org.opencontainers.image.licenses="MIT"


ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# ───────────────────────────────
# System deps
# ───────────────────────────────
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    git wget curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ───────────────────────────────
# Python env
# ───────────────────────────────
RUN python3 -m venv /chai1_venv
ENV PATH="/chai1_venv/bin:$PATH"

RUN pip install --upgrade pip

# ───────────────────────────────
# Torch (CUDA 11.8 stable)
# ───────────────────────────────
RUN pip install torch==2.6.0 torchvision --index-url https://download.pytorch.org/whl/cu118

# ───────────────────────────────
# Chai dependencies
# ───────────────────────────────
RUN pip install \
    chai-lab \
    numpy \
    scipy \
    pandas \
    tqdm \
    biopython \
    safetensors

# ───────────────────────────────
# CRITICAL FIX: force CPU ESM toggle
# ───────────────────────────────
ENV CHAI_FORCE_ESM_CPU=1
ENV CHAI_DOWNLOADS_DIR=/chai_cache

# ───────────────────────────────
# PATCH: robust ESM device handling (NO GPU crash)
# ───────────────────────────────
RUN python3 - <<'PY'
from pathlib import Path

p = Path("/chai1_venv/lib/python3.10/site-packages/chai_lab/data/dataset/embeddings/esm.py")
t = p.read_text()

# 1. force safe esm_device logic
if "CHAI_FORCE_ESM_CPU" not in t:
    t = t.replace(
        "device)",
        "device if not (device.type == 'cuda' and os.environ.get('CHAI_FORCE_ESM_CPU','0')=='1') else torch.device('cpu'))"
    )

# 2. safety patch for token .to(device)
t = t.replace(".to(device)", ".to(device if not (device.type == 'cuda' and os.environ.get('CHAI_FORCE_ESM_CPU','0')=='1') else torch.device('cpu'))")

p.write_text(t)
print("✔ ESM GPU/CPU fallback patched safely")
PY

# ───────────────────────────────
# cache dir
# ───────────────────────────────
RUN mkdir -p /chai_cache

WORKDIR /app