process fetch {
  queue 'datamover'
  container ''

  output:
  path('all.tsv')

  when: { params.databases.mirtrondb?.run }

  script:
  """
  cp ${params.databases.mirtrondb.remote} all.tsv
  """
}

process parse {
  input:
  path(data)

  output:
  path('*.csv')

  script:
  """
  rnac mirtrondb parse $data .
  """
}

workflow mirtrondb {
  emit: data
  main:
  fetch | parse | set { data }
}
