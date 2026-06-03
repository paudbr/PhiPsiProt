FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG HTTP_PROXY
ARG HTTPS_PROXY
ARG http_proxy
ARG https_proxy
ARG NO_PROXY
ARG no_proxy

ENV HTTP_PROXY=${HTTP_PROXY}
ENV HTTPS_PROXY=${HTTPS_PROXY}
ENV http_proxy=${http_proxy}
ENV https_proxy=${https_proxy}
ENV NO_PROXY=${NO_PROXY}
ENV no_proxy=${no_proxy}

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    wget \
    curl \
    vim \
    build-essential \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh

ENV PATH="/root/.local/bin:$PATH"

WORKDIR /opt

RUN git clone https://github.com/RosettaCommons/RFantibody.git

WORKDIR /opt/RFantibody

RUN echo "use_proxy = on" > /root/.wgetrc && \
    echo "http_proxy = ${HTTP_PROXY}" >> /root/.wgetrc && \
    echo "https_proxy = ${HTTPS_PROXY}" >> /root/.wgetrc

RUN timeout 1800 bash include/download_weights.sh

RUN uv sync

ENV PATH="/opt/RFantibody/.venv/bin:$PATH"

WORKDIR /data

ENTRYPOINT ["/bin/bash"]