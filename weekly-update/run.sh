#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'

[ -d work/tmp ] || mkdir -p work/tmp
[ ! -e local.config ] || rm local.config

when=$(date +'%Y-%m-%d')

if [[ -d singularity/bind/r2dt/ ]]; then
  rm -r singularity/bind/r2dt/
fi

mkdir -p singularity/bind/r2dt/data
pushd singularity/bind/r2dt/data
wget -O cms.tar.gz http://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.2/cms.tar.gz
tar xf cms.tar.gz
popd

ln -s weekly-update/update.config local.config

make rust

# Download latest version of nextflow
curl --max-time 10 -s https://get.nextflow.io | bash
res=$?
# Load module as fallback
if test "$res" != "0"; then
  module load nextflow-21.10.6-gcc-9.3.0-tkuemwd
  NF="nextflow"
else
  NF="./nextflow"
fi

#$NF -quiet run -with-report "$when-import.html" -profile test import-data.nf
$NF -quiet run -with-report "$when-analyze.html" -profile test analyze.nf
# $NF -quiet run -with-report "$when-precompute.html" -profile prod precompute.nf
# $NF -quiet run -with-report "$when-search.html" -profile prod search-export.nf
