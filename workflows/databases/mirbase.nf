process mirbase {
  when { params.databases.mirbase.run }

  output:
  path('*.{csv,parquet}')

  """
  scp $params.databases.mirbase.remote mirbase.json
  rnac mirbase parse mirbase.json .
  """
}
