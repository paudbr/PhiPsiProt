FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY pyrosetta-*.whl /tmp/

RUN pip3 install --break-system-packages numpy \
    && pip3 install --break-system-packages /tmp/pyrosetta-*.whl \
    && rm /tmp/pyrosetta-*.whl

WORKDIR /data
CMD  ["/bin/bash"]
