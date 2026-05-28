FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    pymol \
    xvfb \
    python3 \
    && rm -rf /var/lib/apt/lists/*

CMD ["pymol"]
