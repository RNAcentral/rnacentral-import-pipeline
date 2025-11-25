FROM python:3.11.14-trixie

ENV RNA=/rna

WORKDIR $RNA

RUN apt update
RUN apt upgrade -y

# Install all required packages
RUN apt install -y \
    bedtools \
    ca-certificates \
    curl \
    default-mysql-client \
    devscripts \
    freetds-dev \
    gawk \
    gcc \
    git \
    gzip \
    hmmer \
    jq \
    lftp \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libsqlite3-dev \
    libssl-dev \
    libxml2-utils \
    libxml2-dev \
    libzip-dev \
    moreutils \
    mysql-common \
    openssl \
    pandoc \
    patch \
    pgloader \
    postgresql-17 \
    postgresql-client-17 \
    procps \
    python3 \
    python3-dev \
    python3-pip \
    rsync \
    sbcl \
    tabix \
    tar \
    time \
    unzip \
    zlib1g-dev\
    wget


# Install Infernal
RUN \
    cd $RNA/ && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar -xvzf infernal-1.1.2.tar.gz && \
    rm infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && \
    ./configure --prefix=$RNA/infernal-1.1.2 && \
    make && \
    make install && \
    cd easel && \
    make install

# Install blat
RUN \
    wget https://hgwdev.gi.ucsc.edu/~kent/exe/linux/blatSuite.38.zip -O blat.zip && \
    unzip blatSuite.38.zip -d blat_suite && \
    rm blat.zip


# Install seqkit
RUN \
    mkdir seqkit && \
    cd seqkit && \
    wget https://github.com/shenwei356/seqkit/releases/download/v0.10.0/seqkit_linux_amd64.tar.gz && \
    tar xvf seqkit_linux_amd64.tar.gz && \
    rm seqkit_linux_amd64.tar.gz

# Install ribovore
RUN git clone https://github.com/nawrockie/epn-ofile.git && cd epn-ofile && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/epn-options.git && cd epn-options && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/epn-test.git && cd epn-test && git fetch && git fetch --tags && git checkout ribovore-0.40
RUN git clone https://github.com/nawrockie/ribovore.git && cd ribovore && git fetch && git fetch --tags && git checkout ribovore-0.40

# Install htslib
RUN \
    wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar -jxf htslib-1.18.tar.bz2 && \
	rm htslib-1.18.tar.bz2 && \
	cd htslib-1.18 && \
	make && \
    make install

# Install samtools
RUN \
    wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar jxf samtools-1.18.tar.bz2 && \
	rm samtools-1.18.tar.bz2 && \
	cd samtools-1.18 && \
	make && \
    make install

# Install python requirements
ENV RNACENTRAL_IMPORT_PIPELINE="$RNA/rnacentral-import-pipeline"

# Install useful pip version
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
## Add uv install directory to the front of the path
ENV PATH="/root/.local/bin:$PATH"

COPY pyproject.toml $RNACENTRAL_IMPORT_PIPELINE/pyproject.toml
COPY uv.lock $RNACENTRAL_IMPORT_PIPELINE/uv.lock

WORKDIR "$RNA/rnacentral-import-pipeline"
RUN uv sync --no-editable --frozen
ENV PATH="$RNA/rnacentral-import-pipeline/.venv/bin:$PATH"
RUN python3 -m nltk.downloader words

## Download Rust toolchain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

COPY utils ./utils
COPY Makefile Makefile
COPY Cargo.toml Cargo.toml
COPY Cargo.lock Cargo.lock
ENV PATH="$PATH:/root/.cargo/bin"
ENV CARGO_NET_GIT_FETCH_WITH_CLI=true
RUN  make rust

WORKDIR $RNA

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
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE:$PATH"

ENTRYPOINT ["/bin/bash"]
