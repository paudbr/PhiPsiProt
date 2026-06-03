FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-matplotlib \
    procps \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data

CMD ["/bin/bash"]
