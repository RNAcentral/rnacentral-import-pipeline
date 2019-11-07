FROM gcc:8

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
RUN pip3 install --upgrade pip && \
    pip3 install -r $RNACENTRAL_IMPORT_PIPELINE/requirements.txt

RUN python3 -m textblob.download_corpora

# Copy everything that auto-traveler needs into place. This mimicks the same
# directory structure as the auto-traveler Dockerfile.
WORKDIR /
ENV RIBODIR="$RNA/ribovore"
ENV EPNOPTDIR="$RNA/epn-options"
ENV EPNOFILEDIR="$RNA/epn-ofile"
ENV EPNTESTDIR="$RNA/epn-test"
ENV JIFFY_SCRIPTS="$RNA/jiffy-infernal-hmmer-scripts"
ENV TRAVELER="$RNA/traveler"
ENV RNA_STRUCTURE="$RNA/RNAstructure"
ENV RSCAPE="$RNA/rscape"

COPY --from=rnacentral/auto-traveler:dev $EPNOFILEDIR $EPNOFILEDIR
COPY --from=rnacentral/auto-traveler:dev $EPNOPTDIR $EPNOPTDIR
COPY --from=rnacentral/auto-traveler:dev $EPNTESTDIR $EPNTESTDIR
COPY --from=rnacentral/auto-traveler:dev $RIBODIR $RIBODIR
COPY --from=rnacentral/auto-traveler:dev $JIFFY_SCRIPTS $JIFFY_SCRIPTS
COPY --from=rnacentral/auto-traveler:dev $TRAVELER $TRAVELER
COPY --from=rnacentral/auto-traveler:dev $RNA_STRUCTURE $RNA_STRUCTURE
COPY --from=rnacentral/auto-traveler:dev $RSCAPE $RSCAPE

ENV PATH="$JIFFY_SCRIPTS:$PATH"
ENV PATH="$AUTO_TRAVELER_PY:$PATH"
ENV PATH="$RSCAPE/bin:$PATH"
ENV PATH="$RIBODIR:$PATH"
ENV PATH="$RNA_STRUCTURE/exe:$PATH"
ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"

# Install auto-traveler and related data
WORKDIR /
ENV AUTO_TRAVELER_PY="$RNA/auto-traveler"
RUN \
    git clone https://github.com/RNAcentral/auto-traveler.git $AUTO_TRAVELER_PY && \
    cd $AUTO_TRAVELER_PY && \
    git checkout single-entry-point && \
    git pull origin single-entry-point && \
    /usr/local/bin/pip3 install -r requirements.txt

WORKDIR $AUTO_TRAVELER_PY
RUN python3 ./auto-traveler.py setup

ENV PATH="$AUTO_TRAVELER_PY:$PATH"

WORKDIR $RNA

# Setup environmental variables
ENV PERL5LIB="$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:/usr/bin/env:$PERL5LIB"

ENV RIBOTIMEDIR="/usr/bin"
ENV RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV DATAPATH="$RNA_STRUCTURE/data_tables/"

ENV PATH="$TRAVELER/bin:$PATH"
ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"
ENV PATH="$RNA/blatSrc/bin:$PATH"
ENV PATH="$RNA/seqkit:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE:$PATH"
