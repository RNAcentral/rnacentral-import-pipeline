process snodb {
  output:
  path('*.{csv,parquet}')

  script:
  """
  scp $params.databases.snodb.remote snodb.json
  rnac snodb parse snodb.json .
  """
}
