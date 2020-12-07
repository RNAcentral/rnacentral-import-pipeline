process mirgenedb {
  when { params.databases.mirgenedb.run }

  output:
  path('*.csv')

  """
  scp $params.databases.mirgenedb.remote mirgenedb.json
  rnac external mirgenedb mirgenedb.json .
  """
}
