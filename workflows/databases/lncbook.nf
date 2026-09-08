process lncbook {
  output:
  path('*.{csv,parquet}')

  when: params.databases.lncbook?.run

  script:
  """
  wget -O lncbook.json.gz ${params.databases.lncbook.remote}
  gzip -d lncbook.json.gz
  rnac lncbook parse lncbook.json .
  """
}
