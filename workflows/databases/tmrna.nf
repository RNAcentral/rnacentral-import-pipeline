process tmrna {
  when: { params.databases.tmrna.run }

  output:
  path('*.csv')

  """
  cp $params.databases.tmrna.data tmrna.tsv
  rnac tmrna parse tmrna.tsv .
  """
}
