process lncipedia {
  memory '5GB'

  when:
  params.databases.lncipedia.run == true

  output:
  path('*.csv')

  """
  curl ${params.databases.lncipedia.remote} > lncipedia.json
  rnac lncipedia parse lncipedia.json .
  """
}
