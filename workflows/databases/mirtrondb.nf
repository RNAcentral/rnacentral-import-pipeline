process mirtrondb {
  when: { params.databases.mirtrondb.run }

  output:
  path('*.csv')

  """
  cp ${params.databases.mirtrondb.remote} all.tsv
  rnac mirtrondb parse all.tsv .
  """
}
