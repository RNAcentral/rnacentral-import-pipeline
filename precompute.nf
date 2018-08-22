#!/usr/bin/env nextflow

process find_chunks {
  output:
  stdout raw_ranges

  """
  rnac upi-ranges ${params.precompute.max_entries}
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/precompute/query.sql'))
  .set { ranges }

process precompute_range_query {
  beforeScript 'slack db-work precompute-range || true'
  afterScript 'slack db-done precompute-range || true'
  maxForks params.precompute.maxForks

  input:
  set val(min), val(max), file(query) from ranges

  output:
  file 'raw-precompute.json' into precompute_raw

  """
  psql --variable min=$min --variable max=$max -f "$query" '$PGDATABASE' > raw-precompute.json
  """
}

process precompute_range {
  memory '4 GB'

  input:
  file(raw) from precompute_raw

  output:
  file 'precompute.csv' into precompute_results
  file 'qa.csv' into qa_results

  """
  rnac precompute from-file $raw
  """
}

process load_precomputed_data {
  beforeScript 'slack db-work loading-precompute || true'
  afterScript 'slack db-done loading-precompute || true'

  input:
  file('result*.csv') from precompute_results.collect()
  file('qa*.csv') from qa_results.collect()
  file pre_ctl from Channel.fromPath('files/precompute/load.ctl')
  file qa_ctl from Channel.fromPath('files/precompute/qa.ctl')
  file post from Channel.fromPath('files/precompute/post-load.sql')

  // Use cp to copy the ctl files and not stageInMode because we don't want to
  // have to copy a large number of large files before running this. We only
  // need a copy of these two files in the local directory. This is because
  // pgloader will look at the location of the end of symlink when trying to
  // find files and not the current directory or where the symlink lives.
  """
  cp $pre_ctl _pre.ctl
  cp $qa_ctl _qa.ctl
  cat result*.csv | pgloader _pre.ctl
  cat qa*.csv | pgloader _qa.ctl
  psql -f $post "$PGDATABASE"
  """
}
