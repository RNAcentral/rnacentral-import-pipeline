FROM python:3.8-bullseye

ENV RNA=/rna \
    PYTHONFAULTHANDLER=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=random \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=off \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    POETRY_VERSION=1.1.12 \
    POETRY_NO_INTERACTION=1

WORKDIR $RNA

RUN apt-get update
RUN apt-get upgrade -y

# Install all required packages
RUN apt-get install -y \
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
    libsqlite3-dev \
    libssl1.1 \
    libxml2-utils \
    libzip-dev \
    moreutils \
    mysql-common \
    openssl \
    pandoc \
    patch \
    pgloader \
    postgresql-client-11 \
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
    wget https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip && \
    unzip blatSrc35.zip && \
    rm blatSrc35.zip && \
    cd blatSrc && \
    mkdir bin && \
    make MACHTYPE=x86_64 BINDIR=$PWD/bin

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

# Install poetry
RUN pip install "poetry==$POETRY_VERSION"

COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false
RUN poetry install --no-dev --no-root --no-interaction --no-ansi

# Install python requirements
ENV RNACENTRAL_IMPORT_PIPELINE "$RNA/rnacentral-import-pipeline"

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
ENV PATH="$RNA/blatSrc/bin:$PATH"
ENV PATH="$RNA/seqkit:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE:$PATH"

ENTRYPOINT ["/bin/bash"]
