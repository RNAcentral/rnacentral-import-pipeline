process tmrna {
  when: { params.databases.tmrna?.run }

  output:
  path('*.{csv,parquet}')

  script:
  """
  cp $params.databases.tmrna.data tmrna.tsv
  rnac tmrna parse tmrna.tsv .
  """
}
