process hgnc {
  output:
  path('*.{csv,parquet}')

  when: params.databases.hgnc?.run

  script:
  """
  wget -O raw.json $params.databases.hgnc.remote
  rnac hgnc map raw.json
  """
}
