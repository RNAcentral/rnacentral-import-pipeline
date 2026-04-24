process lncipedia {
  memory '5GB'

  when:
  params.databases.lncipedia.run == true

  output:
  path('*.{csv,parquet}')

  """
  wget -O lncipedia.json ${params.databases.lncipedia.remote}
  rnac lncipedia parse lncipedia.json .
  """
}
