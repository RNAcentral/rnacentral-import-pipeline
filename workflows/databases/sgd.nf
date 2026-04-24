process sgd {
  when: { params.databases.sgd.run }

  output:
  path('*.{csv,parquet}')

  """
  wget -O sgd.json ${params.databases.sgd.remote}
  rnac sgd parse sgd.json .
  """
}
