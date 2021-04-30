#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'
export SINGULARITY_TMPDIR='/scratch'

[ -d work/tmp] || mkdir -p work/tmp
[ ! -e local.config ] || rm local.config

when=$(date +'%Y-%m-%d')

ln -s local.config weekly-update/update.config

./nextflow -quiet run -with-report "$when-import.html" -profile prod import-data.nf
./nextflow -quiet run -with-report "$when-search.html" -profile prod search-export.nf
