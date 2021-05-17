#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'
export SINGULARITY_TMPDIR='/scratch'
export PATH="/nfs/software/singularity/3.5.0/bin:$HOME/.cargo/bin:$PATH"

[ -d work/tmp ] || mkdir -p work/tmp
[ ! -e local.config ] || rm local.config

when=$(date +'%Y-%m-%d')

ln -s weekly-update/update.config local.config 

make rust

./nextflow -quiet run -with-report "$when-import.html" -profile prod import-data.nf
./nextflow -quiet run -with-report "$when-precompute.html" -profile prod precompute.nf
./nextflow -quiet run -with-report "$when-search.html" -profile prod search-export.nf
