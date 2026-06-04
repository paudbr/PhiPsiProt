FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    make \
    gcc \
    python3 \
    procps \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/Discngine/fpocket.git /opt/fpocket \
    && cd /opt/fpocket \
    && make \
    && ln -s /opt/fpocket/bin/fpocket /usr/local/bin/fpocket

WORKDIR /data

CMD ["/bin/bash"]
