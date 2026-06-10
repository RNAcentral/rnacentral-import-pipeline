process snodb {
  output:
  path('*.csv')

  script:
  """
  scp $params.databases.snodb.remote snodb.json
  rnac snodb parse snodb.json .
  """
}
