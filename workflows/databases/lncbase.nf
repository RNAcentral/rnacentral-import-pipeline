process lncbase {
  output:
  path('*.{csv,parquet}')

  when: params.databases.lncbase?.run

  script:
  """
  cp ${params.databases.lncbase.remote} lncbase.json
  rnac lncbase parse lncbase.json .
  """
}
