process sgd {
  output:
  path('*.csv')

  when: params.databases.sgd?.run

  script:
  """
  wget -O sgd.json ${params.databases.sgd.remote}
  rnac sgd parse sgd.json .
  """
}
