FROM ubuntu:24.04

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
    python3-dev \
    python3-pip \
    cmake \
    ninja-build \
    build-essential \
    clang \
    llvm \    
    git \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

ENTRYPOINT ["/bin/bash"]
