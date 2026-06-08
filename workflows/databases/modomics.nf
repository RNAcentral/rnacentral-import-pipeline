process fetch {
  queue 'datamover'
  container ''

  output:
  path('modomics.json')

  script:
  """
  cp ${params.databases.modomics.remote} modomics.json
  """
}

process parse {
  input:
  path(data)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac modomics parse $data .
  """
}

workflow modomics {
  emit: data
  main:
  fetch | parse | set { data }
}
