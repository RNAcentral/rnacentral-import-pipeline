process japonicusdb {
  when: { params.databases.japonicusdb.run }

  output:
  path('*.{csv,parquet}')

  """
  wget -O japonicusdb.json ${params.databases.japonicusdb.remote}
  rnac japonicusdb parse japonicusdb.json .
  """
}
