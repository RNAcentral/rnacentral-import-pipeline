process mirbase {
  output:
  path('*.csv')

  script:
  """
  scp $params.databases.mirbase.remote mirbase.json
  rnac mirbase parse mirbase.json .
  """
}
