# Multi-stage build for Infernal 1.1.2
# This container provides the Infernal bioinformatics tool suite

# Stage 1: Builder - compiles Infernal from source
FROM debian:trixie-slim AS builder
ARG INFERNAL_VERSION

ENV RNA=/rna
WORKDIR $RNA

# Install minimal build dependencies
RUN apt update && apt install -y gcc make curl && rm -rf /var/lib/apt/lists/*

# Download and build Infernal
RUN curl -OL http://eddylab.org/infernal/infernal-${INFERNAL_VERSION}.tar.gz && \
    tar -xvzf infernal-${INFERNAL_VERSION}.tar.gz && \
    rm infernal-${INFERNAL_VERSION}.tar.gz && \
    cd infernal-${INFERNAL_VERSION} && \
    ./configure --prefix=$RNA/infernal-${INFERNAL_VERSION} && \
    make && \
    make install && \
    cd easel && \
    make install

# Stage 2: Final - minimal runtime image
FROM debian:trixie-slim
ARG INFERNAL_VERSION
ENV RNA=/rna
WORKDIR $RNA

# Copy built Infernal from builder stage
COPY --from=builder $RNA/infernal-${INFERNAL_VERSION} $RNA/infernal-${INFERNAL_VERSION}

# Set PATH
ENV PATH="$RNA/infernal-${INFERNAL_VERSION}/bin:$PATH"

# Test that Infernal is working
RUN cmscan -h > /dev/null 2>&1

ENTRYPOINT ["/bin/bash"]
