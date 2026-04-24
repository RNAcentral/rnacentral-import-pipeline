process fetch {
  when { params.databases.circatlas.run }
  queue 'datamover'
  container ''

  output:
  path('circatlas.json')

  """
  cp ${params.databases.circatlas.remote} circatlas.json
  """
}

process parse {
  when { params.databases.circatlas.run }

  input:
  path(data)

  output:
  path('*.{csv,parquet}')

  """
  rnac circatlas parse $data .
  """
}

workflow circatlas {
  emit: data
  main:
  fetch | parse | set { data }
}
