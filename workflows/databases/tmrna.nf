process tmrna {
  when: { params.databases.tmrna.run }

  output:
  path('*.{csv,parquet}')

  """
  cp $params.databases.tmrna.data tmrna.tsv
  rnac tmrna parse tmrna.tsv .
  """
}
