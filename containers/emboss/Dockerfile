FROM debian:trixie

RUN apt-get update && \
    apt-get install -y --no-install-recommends emboss procps && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data
