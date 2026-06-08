process japonicusdb {
  when: { params.databases.japonicusdb.run }

  output:
  path('*.{csv,parquet}')

  script:
  """
  wget -O japonicusdb.json ${params.databases.japonicusdb.remote}
  rnac japonicusdb parse japonicusdb.json .
  """
}
