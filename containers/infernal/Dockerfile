# Multi-stage build for Infernal 1.1.2
# This container provides the Infernal bioinformatics tool suite

# Stage 1: Builder - compiles Infernal from source
FROM debian:trixie-slim AS builder

ENV RNA=/rna
WORKDIR $RNA

# Install minimal build dependencies
RUN apt update && apt install -y gcc make curl && rm -rf /var/lib/apt/lists/*

# Download and build Infernal
RUN curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar -xvzf infernal-1.1.2.tar.gz && \
    rm infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && \
    ./configure --prefix=$RNA/infernal-1.1.2 && \
    make && \
    make install && \
    cd easel && \
    make install

# Stage 2: Final - minimal runtime image
FROM debian:trixie-slim

ENV RNA=/rna
WORKDIR $RNA

# Copy built Infernal from builder stage
COPY --from=builder $RNA/infernal-1.1.2 $RNA/infernal-1.1.2

# Set PATH
ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"

# Test that Infernal is working
RUN cmscan -h > /dev/null 2>&1

ENTRYPOINT ["/bin/bash"]
