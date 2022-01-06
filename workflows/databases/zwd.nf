process zwd {
  when: { params.databases.zwd.run }

  input:
  path(context)

  output:
  path('*.csv')

  """
  cp $params.databases.zwd.remote zwd.json
  rnac zwd parse $context zwd.json .
  """
}
