process lncipedia {
  when: { params.databases.lncipedia.run }

  output:
  path('*.csv')

  """
  curl ${params.databases.lncipedia.remote} > lncipedia.json
  rnac external lncipedia lncipedia.json .
  """
}
