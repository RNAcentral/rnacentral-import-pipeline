process lncipedia {
  when: { params.databases.lncipedia.run }
  memory { params.databases.lncipedia.directives.memory }

  output:
  path('*.csv')

  """
  curl ${params.databases.lncipedia.remote} > lncipedia.json
  rnac lncipedia parse lncipedia.json .
  """
}
