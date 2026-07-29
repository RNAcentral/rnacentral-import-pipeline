process tmrna {
  output:
  path('*.csv')

  when: params.databases.tmrna?.run

  script:
  """
  cp $params.databases.tmrna.data tmrna.tsv
  rnac tmrna parse tmrna.tsv .
  """
}
