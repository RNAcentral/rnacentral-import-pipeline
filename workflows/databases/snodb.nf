process snodb {
  when { params.databases.snodb.run }

  output:
  path('*.csv')

  """
  scp $params.databases.snodb.remote snodb.json
  rnac external snodb snodb.json .
  """
}
