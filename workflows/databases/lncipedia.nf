process lncipedia {
  memory '5GB'

  output:
  path('*.csv')

  when: params.databases.lncipedia?.run

  script:
  """
  wget -O lncipedia.json ${params.databases.lncipedia.remote}
  rnac lncipedia parse lncipedia.json .
  """
}
