process lncbase {
  output:
  path('*.csv')

  when: params.databases.lncbase?.run

  script:
  """
  cp ${params.databases.lncbase.remote} lncbase.json
  rnac lncbase parse lncbase.json .
  """
}
