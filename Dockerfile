FROM python:3.7.6-buster

ENV RNA /rna

WORKDIR $RNA

RUN apt-get update
RUN apt-get upgrade -y

# Install all required packages
RUN apt-get install -y \
    gawk \
    bedtools \
    ca-certificates \
    curl \
    default-mysql-client \
    devscripts \
    freetds-dev \
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
    postgresql-11 \
    postgresql-client-11 \
    procps \
    python3 \
    python3-dev \
    python3-pip \
    rsync \
    sbcl \
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

# Install useful pip version
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py

# Install python requirements
ENV RNACENTRAL_IMPORT_PIPELINE "$RNA/rnacentral-import-pipeline"

ADD requirements.txt $RNACENTRAL_IMPORT_PIPELINE/requirements.txt
RUN pip3 install --upgrade pip
RUN pip3 install -r $RNACENTRAL_IMPORT_PIPELINE/requirements.txt

RUN python3 -m textblob.download_corpora

WORKDIR /

WORKDIR $RNA

# Setup environmental variables
ENV PERL5LIB="/usr/bin/env:$PERL5LIB"

ENV RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin"

ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"
ENV PATH="$RNA/blatSrc/bin:$PATH"
ENV PATH="$RNA/seqkit:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE:$PATH"

ENTRYPOINT ["/bin/bash"]
