process japonicusdb {
  output:
  path('*.csv')

  when: { params.databases.japonicusdb?.run }

  script:
  """
  wget -O japonicusdb.json ${params.databases.japonicusdb.remote}
  rnac japonicusdb parse japonicusdb.json .
  """
}
