#!/usr/bin/env nextflow

process find_chunks {
  output:
  stdout raw_ranges

  """
  rnac precompute ranges ${params.precompute.max_entries}
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/precompute/query.sql'))
  .into { ranges }

process precompute_range {
  input:
  set val(min), val(max), file(query) from ranges

  output:
  file 'result.tsv' into precompute_results

  """
  psql --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" |\
  rnac precompute from-file - result.tsv
  """
}

process load_precomputed_data {
  input:
  file('result.tsv*') from precompute_results.collect()
  file ctl from Channel.fromPath('files/precompute/load.ctl')

  """
  cat result.tsv* | split -dC ${params.precompute.load_size} - chunk_

  pgloader $ctl
  """
}
