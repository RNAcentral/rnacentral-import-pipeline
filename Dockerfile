FROM gcc:4.9

ENV RNA /rna
ENV RNACENTRAL_IMPORT_PIPELINE "$RNA/rnacentral-import-pipeline"
RUN mkdir $RNA
WORKDIR $RNA

RUN apt-get install curl ca-certificates
RUN curl https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add -

RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ jessie-pgdg main" > /etc/apt/sources.list.d/pgdg.list

RUN apt-get update
RUN apt-get upgrade -y

# Install all required packages
RUN apt-get install -y \
    bedtools \
    curl \
    gcc \
    git \
    gzip \
    hmmer \
    jq \
    lftp \
    moreutils \
    mysql-client \
    mysql-common \
    openssl \
    pgloader \
    postgresql-9.5 \
    python \
    tar \
    unzip \
    wget

# Install Infernal
RUN \
    cd $RNA/ && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar -xvzf infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && \
    ./configure --prefix=$RNA/infernal-1.1.2 && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd $RNA/ && \
    rm infernal-1.1.2.tar.gz

# Install ribotyper
RUN git clone https://github.com/nawrockie/epn-ofile.git && cd epn-ofile && git checkout c34244b2b9e0719c45d964cc08c147aa353532e8
RUN git clone https://github.com/nawrockie/epn-options.git && cd epn-options && git checkout 7acc13384aedbd5efee9a62fcde71d075072b6a6
RUN git clone https://github.com/nawrockie/epn-test.git && cd epn-test && git checkout f4a8a60153906e61bc458fa734ec7070eadf76f9
RUN git clone https://github.com/nawrockie/ribotyper-v1.git && cd ribotyper-v1 && git checkout 4cd7fe30f402edfa4669383a46d603c60ba6f608

# Install jiffy infernal hmmer scripts
RUN \
    git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git && \
    cd jiffy-infernal-hmmer-scripts && \
    git checkout 45d4937385a6b694eac2d7d538e131b59527ce06

RUN \
    cd jiffy-infernal-hmmer-scripts && \
    echo '#!/usr/bin/env perl' | cat - ali-pfam-sindi2dot-bracket.pl | sponge ali-pfam-sindi2dot-bracket.pl

RUN chmod +x $RNA/jiffy-infernal-hmmer-scripts/ali-pfam-sindi2dot-bracket.pl

# Install traveler
RUN \
    git clone https://github.com/davidhoksza/traveler.git && \
    cd traveler && \
    git checkout 0912ed5daab09bb3c38630efaf3643ea38b02dbe && \
    cd $RNA/traveler/src && \
    make build

# Install auto-traveler.py
RUN git clone https://github.com/RNAcentral/auto-traveler.git && cd auto-traveler && git checkout 5ad1002dc9614e0c0a9c85a0d1a1017ee5027fbe

# Install RNAStructure
RUN \
    wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureSource.tgz && \
    tar -xvzf RNAstructureSource.tgz && \
    rm RNAstructureSource.tgz && \
    cd RNAstructure && \
    make all

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

RUN \
    mkdir -p xsv && \
    cd xsv && \
    wget https://github.com/BurntSushi/xsv/releases/download/0.13.0/xsv-0.13.0-x86_64-unknown-linux-musl.tar.gz && \
    tar xf xsv-0.13.0-x86_64-unknown-linux-musl.tar.gz

# Install useful pip version
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python get-pip.py

# Install python requirements
ADD requirements.txt $RNACENTRAL_IMPORT_PIPELINE/requirements.txt
RUN /usr/local/bin/pip install --upgrade pip && \
    /usr/local/bin/pip install -r $RNACENTRAL_IMPORT_PIPELINE/requirements.txt

RUN python -m textblob.download_corpora

# Setup environmental variables
ENV RIBOINFERNALDIR="$RNA/infernal-1.1.2/bin" RIBOEASELDIR="$RNA/infernal-1.1.2/bin"
ENV RIBODIR="$RNA/ribotyper-v1" EPNOPTDIR="$RNA/epn-options" EPNOFILEDIR="$RNA/epn-ofile" EPNTESTDIR="$RNA/epn-test"
ENV PERL5LIB="$RIBODIR:$EPNOPTDIR:$EPNOFILEDIR:$EPNTESTDIR:/usr/bin/env:$PERL5LIB"

ENV DATAPATH="$RNA/RNAstructure/data_tables/"

ENV PATH="$RNA/traveler/bin:$PATH"
ENV PATH="$RIBODIR:$PATH"
ENV PATH="$RNA/infernal-1.1.2/bin:$PATH"
ENV PATH="$RNA/RNAstructure/exe:$PATH"
ENV PATH="$RNA/blatSrc/bin:$PATH"
ENV PATH="$RNA/seqkit:$PATH"
ENV PATH="$RNA/jiffy-infernal-hmmer-scripts:$PATH"
ENV PATH="$RNA/auto-traveler:$PATH"
ENV PATH="$RNA/xsv:$PATH"
ENV PATH="$RNACENTRAL_IMPORT_PIPELINE:$PATH"
