process mirbase {
  output:
  path('*.csv')

  when: { params.databases.mirbase?.run }

  script:
  """
  scp $params.databases.mirbase.remote mirbase.json
  rnac mirbase parse mirbase.json .
  """
}
