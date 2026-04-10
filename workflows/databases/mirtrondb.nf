process fetch {
  when { params.databases.mirtrondb.run }
  queue 'datamover'
  container ''

  output:
  path('all.tsv')

  """
  cp ${params.databases.mirtrondb.remote} all.tsv
  """
}

process mirtrondb {
  when { params.databases.mirtrondb.run }

  input:
  path(data)

  output:
  path('*.csv')

  """
  rnac mirtrondb parse $data .
  """
}

workflow mirtrondb {
  emit: data
  main:
  fetch | mirtrondb | set { data }
}
