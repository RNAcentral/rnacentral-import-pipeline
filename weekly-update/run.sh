#!/usr/bin/env bash
#BSUB -oo weekly_run.out
#BSUB -eo weekly_run.err
#BSUB -M 4096
#BSUB -cwd /hps/nobackup/agb/rnacentral/weekly-run
#BSUB -J "PDBe weekly import"

set -euo pipefail
IFS=$'\n\t'

export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'

[ -d work/tmp ] || mkdir -p work/tmp
[ ! -e local.config ] || rm local.config

when=$(date +'%Y-%m-%d')


ln -s weekly-update/update.config local.config

make rust

# Download latest version of nextflow
curl --max-time 10 -s https://get.nextflow.io | bash
res=$?
# Load module as fallback
if test "$res" != "0"; then
  echo "Using module nextflow..."
  module load nextflow-21.10.6-gcc-9.3.0-tkuemwd
  NF="nextflow"
else
  echo "Using downloaded nextflow..."
  NF="./nextflow"
fi

# Clean up previous run by nextflow
$NF clean -f

rm .nextflow.log

## Run new DB selection workflow - selects DBs based on file changes from remotes
rm -f db_selection.config && touch db_selection.config

$NF -quiet run -profile pg11prod select_databases.nf --import_selection_remotes=weekly-update/weekly_db_remotes

$NF -quiet run -with-report "$when-setup.html" -profile pg11prod --use_datamover prepare-environment.nf
$NF -quiet run -with-report "$when-import.html" -profile pg11prod import-data.nf
$NF -quiet run -with-report "$when-analyze.html" -profile pg11prod analyze.nf
$NF -quiet run -with-report "$when-precompute.html" -profile pg11prod precompute.nf

$NF -quiet run -profile pg11prod report.nf

$NF -quiet run -with-report "$when-search.html" -profile pg11prod export.nf


# Zip up reports and email them to me
tar -cjf reports.tar.bz2 *.html
rm *.html
mail -a reports.tar.bz2 -s "Weekly workflow completion reports" agreen@ebi.ac.uk < .nextflow.log
