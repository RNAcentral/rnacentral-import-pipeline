# Multi-stage build for Rust utilities
# This container provides all 8 Rust command-line tools used by the pipeline

# Stage 1: Builder - compiles Rust utilities AND maturin Python extension
FROM python:3.11.14-trixie AS builder

# Install Rust toolchain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"
ENV CARGO_NET_GIT_FETCH_WITH_CLI=true

# Install maturin for Python extension building
RUN pip install maturin>=1.10

# Copy Rust source and build
WORKDIR /build
COPY utils ./utils
COPY Cargo.toml Cargo.lock ./

# Build all Rust utilities in release mode
RUN cargo build --release && \
    mkdir /rust-bins && \
    find target/release -maxdepth 1 -type f -executable -exec mv {} /rust-bins/ \;

# Build maturin Python extension wheel
WORKDIR /build/utils/genes-preprocessing
RUN maturin build --release -o /wheels

# Stage 2: Final - minimal runtime image with binaries and wheel
FROM debian:trixie-slim

WORKDIR /rna

# Copy compiled Rust binaries from builder
COPY --from=builder /rust-bins/* /rna/bin/

# Copy maturin wheel from builder
COPY --from=builder /wheels/*.whl /rna/wheels/

# Add binaries to PATH
ENV PATH="/rna/bin:$PATH"

# Test that at least one binary is working
RUN json2fasta --help > /dev/null 2>&1

ENTRYPOINT ["/bin/bash"]
