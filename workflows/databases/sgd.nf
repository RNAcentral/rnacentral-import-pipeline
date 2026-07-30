process sgd {
  output:
  path('*.{csv,parquet}')

  when: params.databases.sgd?.run

  script:
  """
  wget -O sgd.json ${params.databases.sgd.remote}
  rnac sgd parse sgd.json .
  """
}
