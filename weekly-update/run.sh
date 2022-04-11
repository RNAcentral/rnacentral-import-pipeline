#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'
export PATH="/nfs/software/singularity/3.5.0/bin:$HOME/.cargo/bin:$PATH"

[ -d work/tmp ] || mkdir -p work/tmp
[ ! -e local.config ] || rm local.config

when=$(date +'%Y-%m-%d')

if [[ -d singularity/bind/r2dt/ ]]; then
  rm -r singularity/bind/r2dt/
fi

mkdir -p singularity/bind/r2dt/data
pushd singularity/bind/r2dt/data
wget -O cms.tar.gz https://www.dropbox.com/s/3ie8kzb8ol658s0/cms.tar.gz?dl=1
tar xf cms.tar.gz
popd

ln -s weekly-update/update.config local.config 

make rust

./nextflow -quiet run -with-report "$when-import.html" -profile prod import-data.nf
./nextflow -quiet run -with-report "$when-precompute.html" -profile prod precompute.nf
./nextflow -quiet run -with-report "$when-search.html" -profile prod search-export.nf
