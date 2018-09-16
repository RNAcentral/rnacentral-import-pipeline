#!/usr/bin/env nextflow

process fetch_crs {
  input:
  val(remote) from Channel.from(params.crs.path)

  output:
  file('*.tsv.gz') into raw_crs mode flatten

  script:
  """
  fetch "$remote" "*.tsv.gz"
  """
}

process process_crs {
  input:
  file(crs) from raw_crs

  output:
  file('complete_features.csv') into processed_crs

  """
  zcat $crs | rnac extra crs - complete_features.csv
  """
}

process import_crs {
  echo true

  input:
  file('complete_features*.csv') from processed_crs.collect()
  file(ctl) from Channel.fromPath('files/import-metadata/crs-features.ctl')

  """
  cp $ctl crs.ctl
  pgloader --on-error-stop crs.ctl
  """
}
