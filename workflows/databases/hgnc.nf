process hgnc {
  when: { params.databases.hgnc.run }

  output:
  path('*.{csv,parquet}')

  """
  wget -O raw.json $params.databases.hgnc.remote
  rnac hgnc map raw.json
  """
}
