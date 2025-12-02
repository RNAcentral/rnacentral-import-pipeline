# Multi-stage build for Samtools + HTSlib 1.18
# This container provides Samtools, HTSlib tools (tabix, bgzip) and libraries

# Stage 1: Builder - compiles HTSlib and Samtools from source
FROM debian:trixie-slim AS builder

WORKDIR /build

# Install build dependencies
RUN apt update && apt install -y \
    gcc \
    make \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    zlib1g-dev \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Build HTSlib first (Samtools depends on it)
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar -jxf htslib-1.18.tar.bz2 && \
    rm htslib-1.18.tar.bz2 && \
    cd htslib-1.18 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# Build Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar jxf samtools-1.18.tar.bz2 && \
    rm samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# Stage 2: Final - minimal runtime image with libraries
FROM debian:trixie-slim

# Install only runtime dependencies (no -dev packages)
RUN apt update && apt install -y \
    libbz2-1.0 \
    liblzma5 \
    libncurses6 \
    libssl3 \
    zlib1g && \
    rm -rf /var/lib/apt/lists/*

# Copy built binaries and libraries from builder
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/
COPY --from=builder /usr/local/bin/tabix /usr/local/bin/
COPY --from=builder /usr/local/bin/bgzip /usr/local/bin/
COPY --from=builder /usr/local/lib/libhts* /usr/local/lib/

# Run ldconfig to register shared libraries
RUN ldconfig

# Test that samtools and tools are working
RUN samtools --version && \
    tabix --version

ENTRYPOINT ["/bin/bash"]
