process zfin {
  when: { params.databases.zfin.run.toString().toBoolean() }
  errorStrategy 'retry'
  maxRetries 3

  output:
  path('*.csv')

  """
  rnac zfin fetch $params.databases.zfin.remote zfin.json
  rnac zfin parse zfin.json .
  """
}
