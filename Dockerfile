# Multi-stage Dockerfile for RNAcentral Import Pipeline
# This reduces final image size by ~30-40% by separating build and runtime dependencies

# Build arguments for tool versions
ARG INFERNAL_VERSION=1.1.2
ARG SAMTOOLS_VERSION=1.22.1
ARG RUST_VERSION=latest

# Stage 1: Pull pre-built Infernal container
FROM rnacentral/infernal:${INFERNAL_VERSION} AS infernal

# Stage 2: Pull pre-built Samtools/HTSlib container
FROM rnacentral/samtools:${SAMTOOLS_VERSION} AS samtools

# Stage 3: Pull pre-built Rust utilities container
FROM rnacentral/rust-utils:${RUST_VERSION} AS rust-utils

# Stage 4: Python environment builder
FROM python:3.11.14-trixie AS python-builder

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"

# Copy dependency files and install
WORKDIR /app
COPY pyproject.toml .
COPY uv.lock .
RUN uv sync --no-editable --frozen

# Copy and install maturin wheel from rust-utils
COPY --from=rust-utils /rna/wheels/*.whl /tmp/wheels/
RUN uv pip install /tmp/wheels/*.whl

# Download NLTK data
RUN /app/.venv/bin/python3 -m nltk.downloader words

# Stage 5: Final runtime image
FROM python:3.11.14-trixie

ENV RNA=/rna
WORKDIR $RNA

# Install ONLY runtime dependencies (no gcc, no -dev packages)
RUN apt update && apt upgrade -y && \
    apt install -y \
    bedtools \
    ca-certificates \
    curl \
    default-mysql-client \
    gawk \
    git \
    gzip \
    hmmer \
    jq \
    lftp \
    libbz2-1.0 \
    liblzma5 \
    libncurses6 \
    libssl3 \
    libxml2-utils \
    moreutils \
    mysql-common \
    openssl \
    pandoc \
    pgloader \
    postgresql-17 \
    postgresql-client-17 \
    procps \
    python3 \
    python3-pip \
    rsync \
    sbcl \
    tar \
    time \
    unzip \
    wget \
    zlib1g && \
    rm -rf /var/lib/apt/lists/*

# Copy Infernal from tool container
COPY --from=infernal /rna/infernal-1.1.2 $RNA/infernal-1.1.2

# Copy Samtools + HTSlib from tool container
COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/tabix /usr/local/bin/bgzip /usr/local/bin/
COPY --from=samtools /usr/local/lib/libhts* /usr/local/lib/

# Run ldconfig to register shared libraries
RUN ldconfig

# Install blat (pre-compiled)
RUN wget https://hgwdev.gi.ucsc.edu/~kent/exe/linux/blatSuite.38.zip -O blat.zip && \
    unzip blat.zip -d blat_suite && \
    rm blat.zip

# Install seqkit (pre-compiled)
RUN mkdir seqkit && \
    cd seqkit && \
    wget https://github.com/shenwei356/seqkit/releases/download/v2.10.1/seqkit_linux_amd64.tar.gz && \
    tar xvf seqkit_linux_amd64.tar.gz && \
    rm seqkit_linux_amd64.tar.gz

# Clone ribovore (Perl scripts, no compilation)
RUN git clone https://github.com/nawrockie/epn-ofile.git && \
    cd epn-ofile && git checkout ribovore-0.40 && \
    cd .. && \
    git clone https://github.com/nawrockie/epn-options.git && \
    cd epn-options && git checkout ribovore-0.40 && \
    cd .. && \
    git clone https://github.com/nawrockie/epn-test.git && \
    cd epn-test && git checkout ribovore-0.40 && \
    cd .. && \
    git clone https://github.com/nawrockie/ribovore.git && \
    cd ribovore && git checkout ribovore-0.40

# Copy Python environment from builder
ENV RNACENTRAL_IMPORT_PIPELINE="$RNA/rnacentral-import-pipeline"
COPY --from=python-builder /app/.venv $RNACENTRAL_IMPORT_PIPELINE/.venv

# Copy only essential runtime files (exclude build artifacts, tests, Nextflow files)
# Python package - required for CLI
COPY rnacentral_pipeline/ $RNACENTRAL_IMPORT_PIPELINE/rnacentral_pipeline/

# Python/shell scripts - required for various operations (includes old Rust binaries from git)
COPY bin/ $RNACENTRAL_IMPORT_PIPELINE/bin/

# Copy fresh Rust binaries from rust-utils container (overwrites old binaries from git)
COPY --from=rust-utils /rna/bin/* $RNACENTRAL_IMPORT_PIPELINE/bin/

# Package metadata - required for module imports
COPY pyproject.toml setup-env $RNACENTRAL_IMPORT_PIPELINE/


# Setup environmental variables
ENV PERL5LIB="/usr/bin/env:$PERL5LIB"
ENV RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin"
ENV RIBODIR="$RNA/ribovore"
ENV RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV EPNOPTDIR="$RNA/epn-options"
ENV EPNOFILEDIR="$RNA/epn-ofile"
ENV EPNTESTDIR="$RNA/epn-test"
ENV RIBOTIMEDIR="/usr/bin"
ENV BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib"
ENV PERL5LIB="$BIOEASELDIR:$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:$PERL5LIB"

ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"
ENV PATH="$RNA/blat_suite:$PATH"
ENV PATH="$RNA/seqkit:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE/bin:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE/.venv/bin:$PATH"

ENTRYPOINT ["/bin/bash"]
