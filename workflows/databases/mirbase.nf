process mirbase {
  output:
  path('*.{csv,parquet}')

  when: { params.databases.mirbase?.run }

  script:
  """
  scp $params.databases.mirbase.remote mirbase.json
  rnac mirbase parse mirbase.json .
  """
}
