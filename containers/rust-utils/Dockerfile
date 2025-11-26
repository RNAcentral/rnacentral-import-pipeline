# Multi-stage build for Rust utilities
# This container provides all 8 Rust command-line tools used by the pipeline

# Stage 1: Builder - compiles Rust utilities
FROM python:3.11.14-trixie AS builder

# Install Rust toolchain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"
ENV CARGO_NET_GIT_FETCH_WITH_CLI=true

# Copy Rust source and build
WORKDIR /build
COPY utils ./utils
COPY Cargo.toml Cargo.lock Makefile ./

# Build all Rust utilities in release mode
RUN cargo build --release && \
    mkdir /rust-bins && \
    mv target/release/json2fasta \
       target/release/split-ena \
       target/release/expand-urs \
       target/release/precompute \
       target/release/search-export \
       target/release/ftp-export \
       target/release/json2dfasta \
       target/release/bed-expander \
       /rust-bins/

# Stage 2: Final - minimal runtime image
FROM python:3.11.14-trixie

WORKDIR /rna

# Copy compiled Rust binaries from builder
COPY --from=builder /rust-bins/* /rna/bin/

# Add binaries to PATH
ENV PATH="/rna/bin:$PATH"

# Test that at least one binary is working
RUN json2fasta --help > /dev/null 2>&1

ENTRYPOINT ["/bin/bash"]
