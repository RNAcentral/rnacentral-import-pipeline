process lncipedia {
  when: { params.databases.lncipedia.run }
  memory '5GB'

  output:
  path('*.csv')

  """
  curl ${params.databases.lncipedia.remote} > lncipedia.json
  rnac lncipedia parse lncipedia.json .
  """
}
