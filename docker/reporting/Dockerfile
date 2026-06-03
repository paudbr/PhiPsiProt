FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-pandas \
    python3-matplotlib \
    procps \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages multiqc

WORKDIR /data

CMD ["/bin/bash"]
