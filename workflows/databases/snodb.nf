process snodb {
  output:
  path('*.csv')

  when: { params.databases.snodb?.run }

  script:
  """
  scp $params.databases.snodb.remote snodb.json
  rnac snodb parse snodb.json .
  """
}
