FROM gcc:8

ENV RNA /rna

WORKDIR $RNA

RUN apt-get update -m
RUN apt-get upgrade -y

# Install all required packages
RUN apt-get install -y \
    bedtools \
    ca-certificates \
    cpanminus \
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
    libnet-perl \
    libsqlite3-dev \
    libssl1.1 \
    libwww-perl \
    libxml-simple-perl \
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
    python \
    rsync \
    sbcl \
    tar \
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

ENV NCBI="$RNA/ncbi"
WORKDIR $NCBI
RUN curl -s ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz |\
    tar xzf - && \
    cpanm HTML::Entities && \
    edirect/setup.sh

# Install python requirements
ENV RNACENTRAL_IMPORT_PIPELINE "$RNA/rnacentral-import-pipeline"

ADD requirements.txt $RNACENTRAL_IMPORT_PIPELINE/requirements.txt
RUN /usr/local/bin/pip install --upgrade pip && \
    /usr/local/bin/pip install -r $RNACENTRAL_IMPORT_PIPELINE/requirements.txt

RUN python -m textblob.download_corpora

# Copy everything that auto-traveler needs into place. This mimicks the same
# directory structure as the auto-traveler Dockerfile.
WORKDIR /
ENV RIBODIR="$RNA/ribotyper-v1"
ENV EPNOPTDIR="$RNA/epn-options"
ENV EPNOFILEDIR="$RNA/epn-ofile"
ENV EPNTESTDIR="$RNA/epn-test"
ENV JIFFY_SCRIPTS="$RNA/jiffy-infernal-hmmer-scripts"
ENV TRAVELER="$RNA/traveler"
ENV RNA_STRUCTURE="$RNA/RNAstructure"
ENV RSCAPE="$RNA/rscape_v1.2.3"

COPY --from=rnacentral/auto-traveler:rscape-templates $EPNOFILEDIR $EPNOFILEDIR
COPY --from=rnacentral/auto-traveler:rscape-templates $EPNOPTDIR $EPNOPTDIR
COPY --from=rnacentral/auto-traveler:rscape-templates $EPNTESTDIR $EPNTESTDIR
COPY --from=rnacentral/auto-traveler:rscape-templates $RIBODIR $RIBODIR
COPY --from=rnacentral/auto-traveler:rscape-templates $JIFFY_SCRIPTS $JIFFY_SCRIPTS
COPY --from=rnacentral/auto-traveler:rscape-templates $TRAVELER $TRAVELER
COPY --from=rnacentral/auto-traveler:rscape-templates $RNA_STRUCTURE $RNA_STRUCTURE
COPY --from=rnacentral/auto-traveler:rscape-templates $RSCAPE $RSCAPE

ENV PATH="$JIFFY_SCRIPTS:$PATH"
ENV PATH="$AUTO_TRAVELER_PY:$PATH"
ENV PATH="$RSCAPE/bin:$PATH"
ENV PATH="$RIBODIR:$PATH"
ENV PATH="$RNA_STRUCTURE/exe:$PATH"
ENV PATH="$NCBI/edirect:$PATH"

# Install auto-traveler and related data
WORKDIR /
ENV AUTO_TRAVELER_PY="$RNA/auto-traveler"
RUN git clone https://github.com/RNAcentral/auto-traveler.git $AUTO_TRAVELER_PY

WORKDIR $AUTO_TRAVELER_PY
ARG CACHE_DATE=not_a_date
RUN \
    git checkout rscape-templates && \
    /usr/local/bin/pip install -r requirements.txt && \
    ./auto-traveler.py rrna setup && \
    python utils/generate_model_info.py --cm-library data/cms

ENV PATH="$AUTO_TRAVELER_PY:$PATH"

WORKDIR $RNA

# Setup environmental variables
ENV PERL5LIB="$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:/usr/bin/env:$PERL5LIB"

ENV RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV DATAPATH="$RNA_STRUCTURE/data_tables/"

ENV PATH="$TRAVELER/bin:$PATH"
ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"
ENV PATH="$RNA/blatSrc/bin:$PATH"
ENV PATH="$RNA/seqkit:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE:$PATH"
