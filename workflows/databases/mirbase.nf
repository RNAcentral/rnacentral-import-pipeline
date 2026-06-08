process mirbase {
  output:
  path('*.{csv,parquet}')

  script:
  """
  scp $params.databases.mirbase.remote mirbase.json
  rnac mirbase parse mirbase.json .
  """
}
