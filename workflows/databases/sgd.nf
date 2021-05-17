process sgd {
  when: { params.databases.sgd.run }

  output:
  path('*.csv')

  """
  wget -O sgd.json ${params.databases.sgd.remote}
  rnac sgd parse sgd.json .
  """
}
