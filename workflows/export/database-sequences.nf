#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_dbs {
  input:
  path(query)

  output:
  path('dbs.txt')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > dbs.txt
  """
}

process dump_database {
  tag { db }
  maxForks params.export.ftp.sequences.by_database.max_forks

  input:
  tuple val(db), val(full_descr), path(query)

  output:
  tuple val(db), path('raw.json')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" -v db='%${full_descr}%' "\$PGDATABASE" > raw.json
  """
}

// One dump per database, shared by the FTP by-database FASTA and the sequence
// search export, so the (unindexed ilike) query runs once per database.
workflow database_sequences {
  main:
    channel.fromPath('files/ftp-export/sequences/databases.sql') | find_dbs \
    | splitCsv \
    | combine(channel.fromPath('files/ftp-export/sequences/database-specific.sql')) \
    | dump_database
  emit:
    dump_database.out
}
