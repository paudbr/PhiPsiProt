FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Madrid

ARG HTTP_PROXY=http://10.15.206.126:3128
ARG HTTPS_PROXY=http://10.15.206.126:3128
ENV http_proxy=${HTTP_PROXY} https_proxy=${HTTPS_PROXY}
ENV HTTP_PROXY=${HTTP_PROXY} HTTPS_PROXY=${HTTPS_PROXY}
ENV no_proxy=localhost,127.0.0.1
ENV GIT_SSL_NO_VERIFY=true

RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common curl git ca-certificates tzdata && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && apt-get install -y --no-install-recommends \
    python3.12 python3.12-dev python3.12-venv && \
    rm -rf /var/lib/apt/lists/*


# Instalar pip para python3.12
RUN curl -x http://10.15.206.126:3128 \
    https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py && \
    python3.12 /tmp/get-pip.py && \
    rm /tmp/get-pip.py

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1
RUN update-alternatives --install /usr/bin/pip3 pip3 /usr/local/bin/pip3.12 1

# PyTorch 2.4 CUDA 11.8
RUN pip3 install --proxy http://10.15.206.126:3128 \
    torch==2.4.0+cu118 \
    --index-url https://download.pytorch.org/whl/cu118

# Transformers custom de Biohub
RUN git config --global http.proxy http://10.15.206.126:3128 && \
    git config --global https.proxy http://10.15.206.126:3128 && \
    pip3 install --proxy http://10.15.206.126:3128 \
    "transformers @ git+https://github.com/Biohub/transformers.git@main" \
    --no-build-isolation

# Dependencias de ESMFold2
RUN pip3 install --proxy http://10.15.206.126:3128 \
    einops "biotite>=1.0.0" rdkit msgpack-numpy biopython scikit-learn \
    brotli attrs pandas cloudpathlib httpx tenacity zstd \
    ipywidgets py3dmol pydssp boto3 pygtrie accelerate \
    huggingface_hub scipy

# ESMFold2
# Asegura que python3 y pip3 apuntan a 3.12
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1

# ESMFold2 — NO borrar /tmp/esm porque es editable install
RUN git clone https://github.com/Biohub/esm.git /opt/esm && \
    cd /opt/esm && \
    pip3.12 install --proxy http://10.15.206.126:3128 \
    -e . --no-build-isolation --no-deps

# Verificar
RUN python3.12 -c "from esm.models.esmfold2 import ProteinInput; print('ESMFold2 OK')"