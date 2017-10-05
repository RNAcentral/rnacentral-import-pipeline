FROM centos:6.6

RUN yum install -y \
    curl \
    freetds \
    gcc \
    git \
    httpd \
    httpd-devel \
    libaio \
    libzip \
    mysql-devel \
    nc.x86_64 \
    openssl \
    openssl-devel \
    sbcl \
    tar \
    unzip \
    zlib-devel

RUN mkdir /rnacentral
RUN mkdir /rnacentral/local

ENV LOC /rnacentral/local

# Install Python
RUN \
    cd $LOC && \
    curl -OL http://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz && \
    tar -zxvf Python-2.7.11.tgz && \
    cd Python-2.7.11 && \
    PREFIX=$LOC/python-2.7.11/ && \
    export LD_RUN_PATH=$PREFIX/lib && \
    ./configure --prefix=$PREFIX  --enable-shared && \
    make && \
    make install && \
    cd $LOC && \
    rm -Rf Python-2.7.11 && \
    rm Python-2.7.11.tgz

# Install virtualenv
RUN \
    cd $LOC && \
    curl -OL  https://pypi.python.org/packages/source/v/virtualenv/virtualenv-15.0.0.tar.gz && \
    tar -zxvf virtualenv-15.0.0.tar.gz && \
    cd virtualenv-15.0.0 && \
    $LOC/python-2.7.11/bin/python setup.py install && \
    cd $LOC && \
    rm -Rf virtualenv-15.0.0.tar.gz && \
    rm -Rf virtualenv-15.0.0

# Create virtual environment
RUN \
    cd $LOC && \
    mkdir virtualenvs && \
    cd virtualenvs && \
    $LOC/python-2.7.11/bin/virtualenv rnacentral --python=$LOC/python-2.7.11/bin/python

# Install pgloader
RUN \
  rpm --rebuilddb && \
  yum -y install yum-utils && \
  yum -y install rpmdevtools @"Development Tools" && \
  yum -y install sqlite-devel && \
  yum -y install epel-release && \
  yum install -y sbcl.x86_64 --enablerepo=epel

RUN \
  curl -OL http://downloads.sourceforge.net/project/sbcl/sbcl/1.3.6/sbcl-1.3.6-source.tar.bz2 && \
  tar -xvjf sbcl-1.3.6-source.tar.bz2 && \
  cd sbcl-1.3.6 && \
  ./make.sh --with-sb-thread --with-sb-core-compression --prefix=/usr && \
  sh install.sh  && \
  rpm --rebuilddb && \
  yum -y install freetds-devel && \
  rpmdev-setuptree

RUN \
  cd $LOC && \
  curl -OL https://github.com/dimitri/pgloader/archive/v3.4.1.tar.gz && \
  tar -xvzf v3.4.1.tar.gz && \
  rm v3.4.1.tar.gz && \
  mv pgloader-3.4.1 pgloader && \
  cd pgloader && \
  make pgloader

RUN yum -y install gcc

# Install Python requirements
ADD requirements.txt $RNACENTRAL_IMPORT_PIPELINE
RUN \
    source $LOC/virtualenvs/rnacentral/bin/activate && \
    pip install -r $RNACENTRAL_IMPORT_PIPELINE/requirements.txt

# Define container environment variables
ENV PYTHONPATH luigi:$PYTHONPATH
ENV PATH $LOC/pgloader/build/bin:$PATH

COPY entrypoint.sh /entrypoint.sh
RUN chmod 700 entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
