FROM python:3.8-buster

ENV RNA /rna

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
    postgresql-11 \
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

RUN apt-get install -y \
 gconf-service \
 libasound2 \
 libatk1.0-0 \
 libatk-bridge2.0-0 \
 libc6 \
 libcairo2 \
 libcups2 \
 libdbus-1-3 \
 libexpat1 \
 libfontconfig1 \
 libgcc1 \
 libgconf-2-4 \
 libgdk-pixbuf2.0-0 \
 libglib2.0-0 \
 libgtk-3-0 \
 libnspr4 \
 libpango-1.0-0 \
 libpangocairo-1.0-0 \
 libstdc++6 \
 libx11-6 \
 libx11-xcb1 \
 libxcb1 \
 libxcomposite1 \
 libxcursor1 \
 libxdamage1 \
 libxext6 \
 libxfixes3 \
 libxi6 \
 libxrandr2 \
 libxrender1 \
 libxss1 \
 libxtst6 \
 ca-certificates \
 fonts-liberation \
 libappindicator1 \
 libnss3 \
 lsb-release \
 xdg-utils \
 wget \
 libcairo-gobject2 \
 libxinerama1 \
 libgtk2.0-0 \
 libpangoft2-1.0-0 \
 libthai0 \
 libpixman-1-0 \
 libxcb-render0 \
 libharfbuzz0b \
 libdatrie1 \
 libgraphite2-3 \
 libgbm1


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

# Install useful pip version
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py

# Install python requirements
ENV RNACENTRAL_IMPORT_PIPELINE "$RNA/rnacentral-import-pipeline"

ADD requirements.txt $RNACENTRAL_IMPORT_PIPELINE/requirements.txt
RUN pip3 install --upgrade pip
RUN pip3 install -r $RNACENTRAL_IMPORT_PIPELINE/requirements.txt

RUN python3 -m textblob.download_corpora

# Install pyppeteer dependency on chromium
ENV PYPPETEER_HOME=/pyppeteer
RUN mkdir /pyppeteer && pyppeteer-install

WORKDIR /

COPY openssl/openssl.cnf /etc/ssl/

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
