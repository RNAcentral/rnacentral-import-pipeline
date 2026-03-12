FROM python:3.11-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir attrs==25.3.0 stopfree==0.2.4

WORKDIR /data
