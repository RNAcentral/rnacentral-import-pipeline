process zwd {
  when: { params.databases.zwd.run }

  output:
  path('*.csv')

  """
  curl $params.databases.zwd.remote > zwd.json
  rnac external zwd zwd.json .
  """
}
