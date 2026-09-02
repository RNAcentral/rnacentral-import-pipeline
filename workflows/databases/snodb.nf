process snodb {
  output:
  path('*.{csv,parquet}')

  when: params.databases.snodb?.run

  script:
  """
  scp $params.databases.snodb.remote snodb.json
  rnac snodb parse snodb.json .
  """
}
