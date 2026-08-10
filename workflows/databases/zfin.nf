process zfin {
  errorStrategy 'retry'
  maxRetries 3

  output:
  path('*.{csv,parquet}')

  when: params.databases.zfin?.run

  script:
  """
  rnac zfin fetch $params.databases.zfin.remote zfin.json
  rnac zfin parse zfin.json .
  """
}
