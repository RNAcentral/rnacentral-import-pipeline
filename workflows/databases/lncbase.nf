process lncbase {
  when: { params.databases.lncbase.run }

  output:
  path('*.csv')

  """
  cp ${params.databases.lncbase.remote} lncbase.json
  rnac lncbase parse lncbase.json .
  """
}
