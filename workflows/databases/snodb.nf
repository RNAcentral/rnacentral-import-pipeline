process snodb {
  when { params.databases.snodb.run }

  output:
  path('*.csv')

  script:
  """
  scp $params.databases.snodb.remote snodb.json
  rnac snodb parse snodb.json .
  """
}
