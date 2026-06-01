FROM nvidia/cuda:12.6.0-devel-ubuntu24.04

LABEL org.opencontainers.image.title="gnina" \
      org.opencontainers.image.description="Docker image for gnina (CUDA + RDKit + OpenBabel + Torch)" \
      org.opencontainers.image.version="1.3.2"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

WORKDIR /root

# =========================
# 1. System dependencies
# =========================
RUN apt-get update && apt-get install -y \
    build-essential cmake git wget curl swig \
    libboost-all-dev libeigen3-dev libgoogle-glog-dev \
    libprotobuf-dev protobuf-compiler libhdf5-dev \
    libatlas-base-dev libjsoncpp-dev \
    libboost-program-options-dev libboost-system-dev libboost-filesystem-dev \
    python3-dev python3-pip python3-venv python3-pytest \
    librdkit-dev \
    && rm -rf /var/lib/apt/lists/*

# =========================
# 2. Python virtual env
# =========================
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# =========================
# 3. Python deps (IMPORTANT: pytest included!)
# =========================
RUN pip install --no-cache-dir \
    numpy \
    scikit-image \
    pyquaternion \
    google-api-python-client \
    six \
    pytest \
    torch torchvision torchaudio

# =========================
# 4. OpenBabel
# =========================
RUN git clone https://github.com/openbabel/openbabel.git && \
    cd openbabel && mkdir build && cd build && \
    cmake -DWITH_MAEPARSER=OFF \
          -DWITH_COORDGEN=OFF \
          -DPYTHON_BINDINGS=ON \
          -DRUN_SWIG=ON .. && \
    make -j$(nproc) && make install

# =========================
# 5. gnina build
# =========================
ADD https://github.com/gnina/gnina/commits?per_page=1 latest_commit

RUN git clone https://github.com/gnina/gnina.git && \
    cd gnina && mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install