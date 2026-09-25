#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Single point of extraction for the active sequences. Both the FTP FASTA export
// and the HuggingFace parquet export consume this one dump, so the (expensive)
// database read happens once per release. Change how active sequences are pulled
// here and in files/ftp-export/sequences/active.sql, nowhere else.
process dump_active_sequences {

  input:
    path(query)

  output:
    path('active_sequences.json'), emit: json

  script:
  """
  set -euo pipefail
  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > active_sequences.json
  """
}

workflow active_sequences {
  main:
    Channel.fromPath('files/ftp-export/sequences/active.sql') | dump_active_sequences
  emit:
    dump_active_sequences.out.json
}
