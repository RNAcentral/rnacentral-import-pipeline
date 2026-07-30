process psicquic {
  output:
  path('*.{csv,parquet}')

  when: params.databases.psicquic?.run

  script:
  """
  cp $params.databases.psicquic.remote raw.tsv
  rnac psicquic parse raw.tsv .
  """
}
